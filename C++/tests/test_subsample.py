"""CLI subsampling checks using independent V9 and V10 byte fixtures."""
import pathlib
import struct
import sys
import tempfile
import zlib

sys.dont_write_bytecode = True
from test_v10 import fixture, pack, string, run


def legacy_fixture(path, bad_count=None):
    data = bytearray(b'HIC\0' + pack('IQ', 9, 0) + string('test') + bytes(16))
    data += pack('II', 0, 3)
    for name, length in [('All', 1), ('chrA', 80), ('chrB', 70)]:
        data += string(name) + pack('Q', length)
    data += pack('IIII', 2, 10, 20, 0)
    matrices = []
    for a, b, records in [(0, 0, [(0, 0, 1000)]),
                           (1, 1, [(0, 0, 2), (1, 2, 2)]),
                           (1, 2, [(2, 1, 2)]), (2, 2, [(3, 3, 2)])]:
        pos = len(data)
        meta = bytearray(pack('iii', a, b, 2))
        blocks, patches = [], []
        for bp in [10, 20]:
            # Deliberately distinct totals: 8 at fine vs 16 at coarse.
            raw = pack('iii4Bh', len(records), 0, 0, 1, 0, 0, 1, len(records))
            for x, y, n in records:
                count = bad_count if bad_count is not None and a == 1 else n * (bp // 10)
                raw += pack('hhhf', y, 1, x, count)
            block = zlib.compress(raw)
            meta += string('BP') + pack('i4f4i', 0, 0, 0, 0, 0, bp, 8, 1, 1)
            patches.append(len(meta) + 4)
            meta += pack('iQi', 0, 0, len(block))
            blocks.append(block)
        offset = pos + len(meta)
        for patch, block in zip(patches, blocks):
            struct.pack_into('<Q', meta, patch, offset)
            offset += len(block)
        data += meta + b''.join(blocks)
        matrices.append((a, b, pos, len(meta)))
    footer = pack('I', len(matrices))
    for a, b, pos, size in matrices:
        footer += string(f'{a}_{b}') + pack('Qi', pos, size)
    footer += pack('I', 0)
    struct.pack_into('<Q', data, 8, len(data))
    data += pack('Q', len(footer)) + footer + pack('II', 0, 0)
    path.write_bytes(data)


def check(straw, path):
    def command(*args): return [straw, 'subsample', path, *args]
    full = run(command('--fraction', 1))
    rows = [line.split() for line in full.splitlines()]
    assert rows and all(len(row) == 5 and int(row[4]) > 0 for row in rows)
    assert all(row[0].lower() != 'all' and row[2].lower() != 'all' for row in rows)
    assert run(command('--fraction', 0)) == ''
    sampled = run(command('--fraction', .5, '--seed', 27))
    assert sampled == run(command('--fraction', .5, '--seed', 27))
    limits = {tuple(row[:4]): int(row[4]) for row in rows}
    assert all(0 < int(row[4]) <= limits[tuple(row[:4])]
               for row in map(str.split, sampled.splitlines()))
    for args in [('--fraction', -1), ('--fraction', 1.01), ('--fraction', 'nan'),
                 ('--fraction', '0.5junk'), ('--contacts', '-1'), ('--contacts', '1.5'),
                 ('--fraction', .5, '--contacts', 2), ('--seed', 2),
                 ('--fraction', 1, '--resolution', 17), ('--contacts', 100000),
                 ('--fraction', .5, '--seed', '-1'), ('--fraction', .5, '--fraction', .5)]:
        run(command(*args), ok=False)
    return full


def main():
    straw = sys.argv[1]
    with tempfile.TemporaryDirectory() as tmp:
        path = pathlib.Path(tmp) / 'test.hic'
        legacy_fixture(path)
        full = check(straw, path)
        assert full.splitlines() == ['chrA\t0\tchrA\t0\t2', 'chrA\t10\tchrA\t20\t2',
                                     'chrA\t20\tchrB\t10\t2', 'chrB\t30\tchrB\t30\t2']
        base = [straw, 'subsample', path]
        assert run(base + ['--contacts', 8]) == run(base + ['--fraction', .5])
        assert run(base + ['--contacts', 16]) == full
        assert run(base + ['--contacts', 0]) == ''
        assert 'chrA\t40\tchrB\t20\t4' in run(base + ['--fraction', 1, '--resolution', 20])
        # A count of two must produce 0, 1 AND 2, rather than retaining a whole
        # cell or deterministically scaling it. Exercise many fixed seeds.
        outcomes = set()
        for seed in range(32):
            lines = run(base + ['--fraction', .5, '--seed', seed]).splitlines()
            outcomes.add(next((int(line.split()[4]) for line in lines
                               if line.startswith('chrA\t0\tchrA\t0\t')), 0))
        assert outcomes == {0, 1, 2}, outcomes
        for bad in [1.5, -1, float('inf'), float('nan')]:
            legacy_fixture(path, bad)
            run(base + ['--fraction', 1], ok=False)

        fixture(path, values=(2, 2, 2))
        check(straw, path)
        assert run(base + ['--contacts', 3]) == run(base + ['--fraction', .5])
        # Exercise the recursive beta split (n > 64), including its mean and
        # variance. Broad bounds avoid depending on a platform's STL sequence.
        fixture(path, values=(100, 100, 100))
        draws = []
        for seed in range(128):
            rows = run(base + ['--fraction', .2, '--seed', seed]).splitlines()
            counts = [int(row.split()[4]) for row in rows]
            draws.extend(counts + [0] * (3 - len(counts)))
        mean = sum(draws) / len(draws)
        variance = sum((n - mean) ** 2 for n in draws) / len(draws)
        assert abs(mean - 20) < 1 and 10 < variance < 22, (mean, variance)
        fixture(path, values=((1 << 53) + 1, 1, 5))
        assert str((1 << 53) + 1) in run(base + ['--fraction', 1])
        fixture(path, values=((1 << 64) - 3, 1, 1))
        assert str((1 << 64) - 3) in run(base + ['--fraction', 1])
        run(base + ['--fraction', .5])  # Chunked binomial handles full uint64.
        fixture(path, values=((1 << 64) - 1, 1, 1))
        run(base + ['--contacts', 1], ok=False)  # Coarse aggregation overflows.
        fixture(path, score=True, values=(0x3fa00000, 0x40000000, 0x40000000))
        run(base + ['--fraction', 1], ok=False)
        if len(sys.argv) == 4:
            roundtrip(straw, sys.argv[2], sys.argv[3], pathlib.Path(tmp))
        print('Subsampling: V9/V10, cis/trans, coarse totals, seeds, per-count sampling, validation passed')


def roundtrip(straw, pre, v10, tmp):
    chrom = tmp / 'chrom.sizes'
    chrom.write_text('chrA\t10000\nchrB\t7500\n')
    short = tmp / 'input.short'
    text = ('chrA\t100\tchrA\t200\t20\nchrA\t300\tchrA\t400\t10\n'
            'chrA\t500\tchrB\t600\t12\nchrB\t0\tchrB\t100\t8\n')
    short.write_text(text)
    original = tmp / 'original.v9.hic'
    run([pre, '-f', 'short', '-r', '100,200', short, original, chrom])
    converted = tmp / 'converted.v10.hic'
    run([v10, 'convert', original, converted])
    for source in [original, converted]:
        full = run([straw, 'subsample', source, '--fraction', 1])
        assert sorted(full.splitlines()) == sorted(text.splitlines()), (source, full)
        sampled = run([straw, 'subsample', source, '--fraction', .5, '--seed', 42])
        assert sampled == run([straw, 'subsample', source, '--contacts', 25, '--seed', 42])
        short.write_text(sampled)
        for command in [[pre], [v10, 'pre']]:
            rebuilt = tmp / 'rebuilt.hic'
            run(command + ['-f', 'short', '-r', '100,200', short, rebuilt, chrom])
            exported = run([straw, 'subsample', rebuilt, '--fraction', 1])
            assert sorted(exported.splitlines()) == sorted(sampled.splitlines())
    print('V9 and converted V10: subsample -> short -> V9/V10 round trips passed')


if __name__ == '__main__': main()
