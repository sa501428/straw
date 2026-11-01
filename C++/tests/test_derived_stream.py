#!/usr/bin/env python3
"""Derived reads: block boundaries, bounded spills, exact sums, one source pass."""
import collections
import http.server
import math
import pathlib
import random
import struct
import subprocess
import sys
import tempfile
import threading

sys.dont_write_bytecode = True
from test_v10 import pack, string, var, compress, run


def fixture(path, cells, length, factor, block_size, score=False, trans=False):
    """Independent multi-block writer with deliberately unaligned bin sizes."""
    columns = math.ceil(length / block_size)
    groups = collections.defaultdict(list)
    for x, y, value in cells:
        if trans:
            number = (y // block_size) * columns + x // block_size
        else:
            depth = math.floor(math.log2(1 + (y-x) / (math.sqrt(2) * block_size)))
            number = depth * columns + (x+y) // (2*block_size)
        groups[number].append((x, y, value))
    variable = string('stream-test') + pack('II', 0, 2)
    variable += string('chrA') + pack('Q', length) + string('chrB') + pack('Q', length)
    variable += pack('I', 2) + pack('IBBHI', 1, 0, 1, 0, 0xffffffff)
    variable += pack('IBBHI', factor, 1, 1, 0, 0) + pack('II', 0, 0)
    data = bytearray(88) + variable
    matrix_pos = len(data)
    data += b'H10M' + pack('5I', 1, 0, int(trans), 2, 0) + bytes(152)
    entries = []
    for number, group in sorted(groups.items()):
        group.sort(key=lambda r: (r[1], r[0]))
        x0, y0 = min(r[0] for r in group), min(r[1] for r in group)
        width = max(r[0] for r in group) - x0 + 1
        height = max(r[1] for r in group) - y0 + 1
        positions = [(y-y0)*width + x-x0 for x, y, _ in group]
        ps = var(positions[0]) + b''.join(var(b-a) for a, b in zip(positions, positions[1:]))
        vs = b''.join(pack('f', v) if score else var(v) for _, _, v in group)
        raw = pack('5B3x4IQ2I', 1, 0, 2, int(score), 0, x0, y0, width, height,
                   len(group), len(ps), len(vs)) + ps + vs
        stored = b'H10B' + pack('BBHII', 1, 1, 0, len(raw), number) + compress(raw)
        entries.append((number, len(stored), len(data)))
        data += stored
    index_pos = len(data)
    index = b'H10I' + pack('IQII', 2, 24 + 16*len(entries), len(entries), 0)
    index += b''.join(pack('IIQ', *entry) for entry in entries)
    data += index
    total = sum(r[2] for r in cells)
    sum_bits = int.from_bytes(pack('d', total), 'little') if score else total
    occupied = len({(x//factor, y//factor) for x, y, _ in cells})

    def descriptor(derived):
        bin_size = factor if derived else 1
        return pack('4B3IB3xQQII2IQQ2I', 0, int(derived), 1, int(score), int(derived), bin_size,
                    0 if derived else 0xffffffff, int(not trans), sum_bits,
                    occupied if derived else len(cells), 0x7fc00000, 0x7fc00000,
                    block_size, math.ceil(math.ceil(length/bin_size) / block_size),
                    0 if derived else index_pos, 0 if derived else len(index),
                    0 if derived else len(entries), 0)

    data[matrix_pos+24:matrix_pos+176] = descriptor(False) + descriptor(True)
    footer_pos = len(data)
    data += b'H10F' + pack('IQII', 1, 48, 1, 0) + pack('IIQQ', 0, int(trans), matrix_pos, 176)
    data[:88] = b'HIC\0' + pack('I10Q', 10, 88+len(variable), footer_pos, 48, 0, 0, 0, 0, 0, 0, 0)
    path.write_bytes(data)
    return entries


def expected(cells, factor, score, trans=False, bounds=None, reverse=False):
    sums = collections.defaultdict(float if score else int)
    for x, y, value in sorted(cells, key=lambda r: (r[1], r[0])):
        x, y = x//factor, y//factor
        if reverse:
            x, y = y, x
        if bounds:
            a, b, c, d = bounds
            if not (a <= x < b and c <= y < d) and not (
                    not trans and a <= y < b and c <= x < d):
                continue
        sums[x, y] += value
    return {k: ('s', struct.unpack('<I', pack('f', v))[0]) if score else ('c', v)
            for k, v in sums.items()}


def check(probe, path, factor, cells, score=False, trans=False, bounds=None, reverse=False):
    names = ['chrA', 'chrB' if trans else 'chrA']
    if reverse:
        names.reverse()
    args = [probe, path, 'raw', *names, factor, 'BP']
    if bounds:
        args += bounds
    rows = [r.split() for r in run(args).splitlines()]
    actual = {(int(x), int(y)): (kind, int(value)) for x, y, kind, value in rows}
    assert len(actual) == len(rows), 'a target cell was emitted more than once'
    wanted = expected(cells, factor, score, trans, bounds, reverse)
    assert actual == wanted, (path, factor, score, bounds, reverse,
                              list((actual.items() - wanted.items()))[:5],
                              list((wanted.items() - actual.items()))[:5])


def http_check(probe, path, factor, entries, stream_before=None):
    """Count actual source reads; optionally demand output before a later read."""
    data = path.read_bytes()
    requests = collections.Counter()
    first = threading.Event()
    failures = []

    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self, *_):
            pass

        def do_GET(self):
            start, end = map(int, self.headers['Range'][6:].split('-'))
            requests[start, end-start+1] += 1
            if stream_before == start and not first.wait(5):
                failures.append('no aggregate emitted before reading the final source block')
            self.send_response(206)
            self.send_header('Content-Range', f'bytes {start}-{end}/{len(data)}')
            self.send_header('Content-Length', str(end-start+1))
            self.end_headers()
            self.wfile.write(data[start:end+1])

    server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        url = f'http://127.0.0.1:{server.server_port}/fixture.hic'
        proc = subprocess.Popen([str(probe), url, 'raw-progress', 'chrA', 'chrA', str(factor)],
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            assert proc.stdout.readline() == 'first\n'
            first.set()
            output, error = proc.communicate(timeout=20)
            assert proc.returncode == 0, error
            assert output.startswith('records '), output
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.communicate()
        assert not failures, failures
        for _, size, pos in entries:
            assert requests[pos, size] == 1, (pos, size, requests[pos, size])
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


def main():
    probe = sys.argv[1]
    with tempfile.TemporaryDirectory() as directory:
        path = pathlib.Path(directory) / 'derived.hic'
        rng = random.Random(42)
        # Source and target boundaries cross both rotated cis bands and
        # rectangular trans blocks; last chromosome bins are partial.
        for trans in (False, True):
            for score in (False, True):
                values = [1., -1., 2.**60, -2.**60, 3., 0., -0.] if score else [1, 3, 2**53+7]
                cells = [(x, y, rng.choice(values)) for y in range(83) for x in range(83)
                         if (trans or x <= y) and rng.random() < .3]
                for factor in (2, 5, 11, 97):
                    fixture(path, cells, 83, factor, 7, score, trans)
                    check(probe, path, factor, cells, score, trans)
                    check(probe, path, factor, cells, score, trans, (0, 7, 3, 15))
                    check(probe, path, factor, cells, score, trans, (3, 15, 0, 7), trans)

        # A huge single source block should stream without recursively decoding
        # it for each subregion. This matches the old 65,536-cell regression.
        cells = [(2*i, 2*i+2, 1) for i in range(70000)]
        entries = fixture(path, cells, 140002, 2, 140002)
        check(probe, path, 2, cells)
        if '--http' in sys.argv:
            http_check(probe, path, 2, entries)

        # Many unfinished targets on a single source row force several sorted
        # spill runs. Repeated source cells for each target test merging sums.
        cells = [(x, 280003, 1) for x in range(280000)]
        entries = fixture(path, cells, 280005, 2, 280005)
        check(probe, path, 2, cells)
        if '--http' in sys.argv:
            http_check(probe, path, 2, entries)

        # One target score cell spanning many blocks exceeds the record budget.
        # Cancellation makes summing partial block scores observably incorrect.
        cells = [(x, y, [2.**60, 1., -2.**60, 3.][(x+3*y) % 4])
                 for y in range(300) for x in range(1000)]
        fixture(path, cells, 1201, 1000, 97, score=True, trans=True)
        check(probe, path, 1000, cells, score=True, trans=True)

        # A completed early region must reach the callback before the final
        # source block is even fetched, rather than collecting a whole matrix.
        cells = [(0, 0, 1), (0, 4, 2), (70, 70, 3)]
        entries = fixture(path, cells, 83, 5, 7)
        check(probe, path, 5, cells)
        if '--http' in sys.argv:
            http_check(probe, path, 5, entries, stream_before=entries[-1][2])
    print('Derived streaming: boundaries, exact counts, score order, bounded spills passed')


if __name__ == '__main__':
    main()
