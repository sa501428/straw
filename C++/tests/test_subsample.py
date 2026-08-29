"""CLI subsampling checks using independent V9 and V10 byte fixtures."""
import pathlib
import gzip
import os
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


def read_hbs(path):
    """Independent HBS decoder, including exact record widths."""
    data = gzip.decompress(path.read_bytes())
    assert data[:8] == b'HICBS\0\r\n'
    version, flags, resolution, nchr = struct.unpack_from('<HHII', data, 8)
    assert version == 1 and flags == 0
    offset, chroms = 20, []
    for _ in range(nchr):
        length, = struct.unpack_from('<H', data, offset)
        offset += 2
        name = data[offset:offset+length].decode()
        offset += length
        bp, = struct.unpack_from('<Q', data, offset)
        offset += 8
        chroms.append((name, bp))
    rows = []
    while offset < len(data):
        a, x, b, y, count = struct.unpack_from('<HIHIH', data, offset)
        offset += 14
        if count == 65535:
            count, = struct.unpack_from('<Q', data, offset)
            offset += 8
            assert count >= 65535
        assert count > 0 and a < nchr and b < nchr
        rows.append(f'{chroms[a][0]}\t{x*resolution}\t{chroms[b][0]}\t{y*resolution}\t{count}')
    assert offset == len(data)
    return resolution, chroms, rows


def check_hbs(straw, path):
    output = path.with_name('output.hbs.gz')
    base = [straw, 'subsample', path]
    text = run(base + ['--fraction', .5, '--seed', 42])
    assert run(base + ['--fraction', .5, '--seed', 42, '--output', output]) == ''
    resolution, chroms, rows = read_hbs(output)
    assert resolution == 10 and chroms == [('chrA', 80), ('chrB', 70)]
    assert rows == text.splitlines()
    for args in [['--fraction', 0], ['--contacts', 0]]:
        run(base + args + ['-o', output])
        assert read_hbs(output)[::2] == (10, [])
    dump = [straw, 'dump', 'observed', 'NONE', path, 'BP', 10, output]
    run(dump)
    full = run(base + ['--fraction', 1]).splitlines()
    assert read_hbs(output)[2] == full
    run(dump + [1]) # Existing dump syntax remains accepted.
    assert read_hbs(output)[2] == full
    for filter, cis in [('-inter', False), ('-intra', True)]:
        run(dump + [filter])
        assert read_hbs(output)[2] == [r for r in full if (r.split()[0] == r.split()[2]) == cis]
    for bad in [[straw, 'dump', 'observed', 'VC', path, 'BP', 10, output],
                [straw, 'dump', 'oe', 'NONE', path, 'BP', 10, output],
                [straw, 'dump', 'observed', 'NONE', path, 'FRAG', 1, output]]:
        run(bad, ok=False)
    alias = path.with_name('alias.hbs.gz')
    os.link(path, alias)
    original = path.read_bytes()
    run(base + ['--fraction', 1, '-o', alias], ok=False)
    assert path.read_bytes() == original
    alias.unlink()


def main():
    straw = sys.argv[1]
    with tempfile.TemporaryDirectory() as tmp:
        path = pathlib.Path(tmp) / 'test.hic'
        legacy_fixture(path)
        # All reader versions return contact axes in request order. Exercise
        # both a trans request opposite to file chromosome order and a cis
        # request below the stored upper triangle.
        assert run([straw, 'observed', 'NONE', path, 'chrB', 'chrA', 'BP', 10]) == '10\t20\t2\n'
        assert run([straw, 'observed', 'NONE', path,
                    'chrA:20:20', 'chrA:10:10', 'BP', 10]) == '20\t10\t2\n'
        matrix = run([straw, 'observed', 'NONE', path, 'chrB', 'chrA', 'MATRIX', 10])
        matrix = [list(map(float, row.split())) for row in matrix.splitlines()]
        assert matrix[1][2] == 2
        full = check(straw, path)
        check_hbs(straw, path)
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
        check_hbs(straw, path)
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
        output = path.with_name('exact.hbs.gz')
        fixture(path, values=(65534, 65535, (1 << 53) + 1))
        run(base + ['--fraction', 1, '-o', output])
        assert [int(row.split()[4]) for row in read_hbs(output)[2]] == [65534, 65535, (1 << 53) + 1]
        fixture(path, values=((1 << 64) - 3, 1, 1))
        assert str((1 << 64) - 3) in run(base + ['--fraction', 1])
        run(base + ['--fraction', .5])  # Chunked binomial handles full uint64.
        fixture(path, values=((1 << 64) - 1, 1, 1))
        run(base + ['--contacts', 1], ok=False)  # Coarse aggregation overflows.
        fixture(path, score=True, values=(0x3fa00000, 0x40000000, 0x40000000))
        run(base + ['--fraction', 1], ok=False)
        fixture(path, score=True, values=(0x3f800000, 0x3fa00000, 0x40000000))
        original = output.read_bytes()
        run(base + ['--fraction', 1, '-o', output], ok=False)
        assert output.read_bytes() == original and not list(output.parent.glob('*.tmp-*'))
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
        binary = tmp / 'sampled.hbs.gz'
        run([straw, 'subsample', source, '--fraction', .5, '--seed', 42, '-o', binary])
        assert sorted(read_hbs(binary)[2]) == sorted(sampled.splitlines())
        # Reverse genome order to test name mapping through the HBS table.
        reverse = tmp / 'reverse.sizes'
        reverse.write_text('chrB\t7500\nchrA\t10000\n')
        for command in [[pre], [v10, 'pre']]:
            for genome in [chrom, reverse]:
                rebuilt = tmp / 'rebuilt.hic'
                run(command + ['-r', '100,200', binary, rebuilt, genome])
                exported = run([straw, 'subsample', rebuilt, '--fraction', 1])
                def canonical(text):
                    result = []
                    for line in text.splitlines():
                        a, x, b, y, n = line.split()
                        if a > b: a, x, b, y = b, y, a, x
                        result.append((a, x, b, y, n))
                    return sorted(result)
                assert canonical(exported) == canonical(sampled)
            run(command + ['-r', 50, binary, rebuilt, chrom], ok=False)
    # Exact escaped counts must survive V10; V9 rounds only at its float path.
    source = tmp / 'exact.v10.hic'
    fixture(source, values=(65534, 65535, (1 << 53) + 1))
    chrom.write_text('chrA\t80\nchrB\t70\n')
    binary = tmp / 'exact.hbs.gz'
    run([straw, 'dump', 'observed', 'NONE', source, 'BP', 10, binary])
    for command in [[pre], [v10, 'pre']]:
        rebuilt = tmp / 'exact-rebuilt.hic'
        run(command + ['-f', 'hbs', '-r', '10,20', binary, rebuilt, chrom])
        exported = run([straw, 'subsample', rebuilt, '--fraction', 1])
        counts = [int(line.split()[4]) for line in exported.splitlines()]
        assert counts == [65534, 65535, (1 << 53) + (1 if command[0] == v10 else 0)], counts
    print('V9 and converted V10: subsample -> short -> V9/V10 round trips passed')


if __name__ == '__main__': main()
