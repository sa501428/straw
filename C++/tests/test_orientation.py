#!/usr/bin/env python3
"""Regression tests for the public user-axis contact orientation contract."""

import collections
import pathlib
import subprocess
import sys


def query(straw, fixture, first, second, unit="BP", resolution=2_500_000):
    process = subprocess.run(
        [straw, "observed", "NONE", fixture, first, second, unit, str(resolution)],
        check=True,
        capture_output=True,
        text=True,
    )
    records = []
    for line in process.stdout.splitlines():
        x, y, value = line.split("\t")
        records.append((int(x), int(y), float(value)))
    return collections.Counter(records)


def transpose(records):
    return collections.Counter((y, x, value) for (x, y, value), count in records.items()
                               for _ in range(count))


def main():
    straw = sys.argv[1]
    fixture = str(pathlib.Path(sys.argv[2]).resolve())

    forward = query(straw, fixture, "1", "2")
    reverse = query(straw, fixture, "2", "1")
    assert forward, "legacy fixture did not return interchromosomal contacts"
    assert reverse == transpose(forward), "reversed chromosome query was not transposed"

    upper = query(straw, fixture, "1:0:5000000", "1:10000000:15000000")
    lower = query(straw, fixture, "1:10000000:15000000", "1:0:5000000")
    assert upper, "legacy fixture did not return asymmetric cis contacts"
    assert lower == transpose(upper), "reflected cis query was not transposed"


if __name__ == "__main__":
    main()
