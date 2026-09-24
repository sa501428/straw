#!/usr/bin/env python3
"""Expand binned HBS counts into independently jittered 1-bp contacts."""

from __future__ import annotations

import argparse
import gzip
import os
import random
import struct
import sys
import tempfile
from pathlib import Path
from typing import BinaryIO


MAGIC = b"HICBS\x00\r\n"
HEADER = struct.Struct("<8sHHII")
RECORD = struct.Struct("<HIHIH")
U16 = struct.Struct("<H")
U64 = struct.Struct("<Q")
MAX_HEADER_BYTES = 16 * 1024 * 1024
MAX_U32 = (1 << 32) - 1
FLUSH_BYTES = 4 * 1024 * 1024


def read_exact(stream: BinaryIO, size: int, description: str) -> bytes:
    chunks: list[bytes] = []
    remaining = size
    while remaining:
        chunk = stream.read(remaining)
        if not chunk:
            raise ValueError(f"truncated HBS while reading {description}")
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def read_record_or_eof(stream: BinaryIO) -> bytes | None:
    first = stream.read(1)
    if not first:
        return None
    return first + read_exact(stream, RECORD.size - 1, "contact record")


def read_header(stream: BinaryIO) -> tuple[int, list[tuple[bytes, int]], bytes]:
    fixed = read_exact(stream, HEADER.size, "header")
    magic, version, flags, resolution, chromosome_count = HEADER.unpack(fixed)
    if magic != MAGIC:
        raise ValueError("invalid HBS magic")
    if version != 1:
        raise ValueError(f"unsupported HBS version: {version}")
    if flags != 0:
        raise ValueError(f"unsupported HBS flags: {flags}")
    if not 1 <= resolution <= 0x7FFFFFFF:
        raise ValueError(f"invalid HBS resolution: {resolution}")
    if chromosome_count > 65536:
        raise ValueError(f"too many chromosomes: {chromosome_count}")

    encoded = bytearray(fixed)
    chromosomes: list[tuple[bytes, int]] = []
    names: set[bytes] = set()
    for chromosome_id in range(chromosome_count):
        name_size_bytes = read_exact(stream, U16.size, "chromosome name length")
        name_size = U16.unpack(name_size_bytes)[0]
        if not 1 <= name_size <= 4096:
            raise ValueError(f"invalid chromosome name length for ID {chromosome_id}")
        name = read_exact(stream, name_size, "chromosome name")
        length_bytes = read_exact(stream, U64.size, "chromosome length")
        length = U64.unpack(length_bytes)[0]
        if b"\x00" in name or name in names or length == 0:
            raise ValueError(f"invalid chromosome table entry for ID {chromosome_id}")
        names.add(name)
        chromosomes.append((name, length))
        encoded.extend(name_size_bytes)
        encoded.extend(name)
        encoded.extend(length_bytes)
        if len(encoded) > MAX_HEADER_BYTES:
            raise ValueError("HBS header exceeds 16 MiB")

    return resolution, chromosomes, bytes(encoded)


def offset_position(anchor: int, length: int, max_offset: int, rng: random.Random) -> int:
    if max_offset:
        magnitude = rng.randint(1, max_offset)
        anchor += -magnitude if rng.getrandbits(1) else magnitude
    return min(max(anchor, 0), length - 1)


def jitter_hbs(
    input_path: Path,
    output_path: Path,
    max_offset: int,
    seed: int | None,
    force: bool,
) -> tuple[int, int, int]:
    if input_path.resolve() == output_path.resolve():
        raise ValueError("input and output paths must differ")
    if not input_path.is_file():
        raise ValueError(f"input file does not exist: {input_path}")
    if output_path.exists() and not force:
        raise ValueError(f"output already exists (use --force to replace it): {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rng = random.Random(seed)
    input_records = 0
    output_records = 0
    total_count = 0
    temporary_name: str | None = None

    try:
        with gzip.open(input_path, "rb") as source:
            source_resolution, chromosomes, original_header = read_header(source)
            # Preserve the chromosome table, but advertise one-base-pair bins.
            output_header = bytearray(original_header)
            output_header[12:16] = struct.pack("<I", 1)

            temporary = tempfile.NamedTemporaryFile(
                mode="wb",
                prefix=f".{output_path.name}.",
                suffix=".tmp",
                dir=output_path.parent,
                delete=False,
            )
            temporary_name = temporary.name
            with temporary as raw_output:
                with gzip.GzipFile(
                    filename="",
                    mode="wb",
                    compresslevel=6,
                    fileobj=raw_output,
                    mtime=0,
                ) as destination:
                    destination.write(output_header)
                    buffer = bytearray()

                    while True:
                        packed = read_record_or_eof(source)
                        if packed is None:
                            break
                        chromosome_1, bin_1, chromosome_2, bin_2, count = RECORD.unpack(packed)
                        if chromosome_1 >= len(chromosomes) or chromosome_2 >= len(chromosomes):
                            raise ValueError("contact record contains a chromosome ID outside the table")
                        if count == 65535:
                            count = U64.unpack(read_exact(source, U64.size, "escaped count"))[0]
                            if count < 65535:
                                raise ValueError("noncanonical escaped HBS count")

                        length_1 = chromosomes[chromosome_1][1]
                        length_2 = chromosomes[chromosome_2][1]
                        anchor_1 = bin_1 * source_resolution
                        anchor_2 = bin_2 * source_resolution
                        if anchor_1 > length_1 or anchor_2 > length_2:
                            raise ValueError("contact bin start lies outside its chromosome")
                        # HBS permits a terminal endpoint equal to chromosome length.
                        anchor_1 = min(anchor_1, length_1 - 1)
                        anchor_2 = min(anchor_2, length_2 - 1)

                        input_records += 1
                        total_count += count
                        for _ in range(count):
                            position_1 = offset_position(anchor_1, length_1, max_offset, rng)
                            position_2 = offset_position(anchor_2, length_2, max_offset, rng)
                            if chromosome_1 == chromosome_2 and position_1 > position_2:
                                position_1, position_2 = position_2, position_1
                            if position_1 > MAX_U32 or position_2 > MAX_U32:
                                raise ValueError("1-bp output position exceeds the HBS uint32 bin limit")
                            buffer.extend(
                                RECORD.pack(
                                    chromosome_1,
                                    position_1,
                                    chromosome_2,
                                    position_2,
                                    1,
                                )
                            )
                            output_records += 1
                            if len(buffer) >= FLUSH_BYTES:
                                destination.write(buffer)
                                buffer.clear()
                            if output_records % 10_000_000 == 0:
                                print(
                                    f"Wrote {output_records:,} pseudo-contacts...",
                                    file=sys.stderr,
                                    flush=True,
                                )

                    if buffer:
                        destination.write(buffer)
                raw_output.flush()
                os.fsync(raw_output.fileno())

        os.replace(temporary_name, output_path)
        temporary_name = None
        return source_resolution, input_records, total_count
    finally:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Expand every HBS count into count-1 contacts, independently jitter both "
            "endpoints, and write an HBS file with 1-bp resolution."
        )
    )
    parser.add_argument("input", type=Path, help="input .hbs.gz file")
    parser.add_argument("output", type=Path, help="output .hbs.gz file")
    parser.add_argument(
        "--max-offset",
        type=int,
        default=500,
        help="largest absolute offset in bp; zero disables jitter (default: 500)",
    )
    parser.add_argument("--seed", type=int, help="random seed for reproducible output")
    parser.add_argument("--force", action="store_true", help="replace an existing output file")
    args = parser.parse_args()
    if args.max_offset < 0:
        parser.error("--max-offset must be nonnegative")
    if not str(args.input).endswith(".hbs.gz") or not str(args.output).endswith(".hbs.gz"):
        parser.error("input and output names must end in .hbs.gz")
    return args


def main() -> int:
    args = parse_args()
    try:
        resolution, records, count = jitter_hbs(
            args.input,
            args.output,
            args.max_offset,
            args.seed,
            args.force,
        )
    except (OSError, EOFError, ValueError, struct.error) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1

    print(f"Input resolution: {resolution:,} bp")
    print(f"Input HBS records: {records:,}")
    print(f"Output count-1 records: {count:,}")
    print(f"Created: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
