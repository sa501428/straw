"""Independent fixtures exercising subdivision of derived query buffers.

Regenerate Rust's checked-in copies with:
python3 review_fixture.py ../../../straw-rust/tests/data
"""
import pathlib
import sys
sys.dont_write_bytecode = True
from test_v10 import pack, string, var, compress


def fixture(path, score=False, fanin=False):
    n = 70000  # exceeds the 65,536-entry aggregation buffer
    length, target = (1000, 1000) if fanin else (2*n+2, 2)
    width, height = (300, 234) if fanin else (2*n+2, 2*n+2)
    positions = list(range(n)) if fanin else [2*i + (2*i+2)*width for i in range(n)]
    ps = var(positions[0]) + b''.join(var(b-a) for a, b in zip(positions, positions[1:]))
    vs = pack('f', 1.) if score else var(1)
    raw = pack('5B3x4IQ2I', 1, 0, 0, int(score), 0, 0, 0, width, height, n, len(ps), len(vs)) + ps + vs
    variable = string('review') + pack('II', 0, 2)
    variable += string('chrA') + pack('Q', length) + string('chrB') + pack('Q', length)
    variable += pack('I', 2) + pack('IBBHI', 1, 0, 1, 0, 0xffffffff) + pack('IBBHI', target, 1, 1, 0, 0)
    variable += pack('II', 0, 0)
    data = bytearray(88) + variable
    mp = len(data)
    data += b'H10M' + pack('5I', 1, 0, int(fanin), 2, 0) + bytes(152)
    bp = len(data)
    stored = b'H10B' + pack('BBHII', 1, 1, 0, len(raw), 0) + compress(raw)
    data += stored
    ip = len(data)
    data += b'H10I' + pack('IQII', 2, 40, 1, 0) + pack('IIQ', 0, len(stored), bp)
    sum_bits = int.from_bytes(pack('d', float(n)), 'little') if score else n
    def descriptor(derived):
        return pack('4B3IB3xQQII2IQQ2I', 0, int(derived), 1, int(score), int(derived), target if derived else 1,
                    0 if derived else 0xffffffff, int(not fanin), sum_bits, 1 if fanin and derived else n,
                    0x7fc00000, 0x7fc00000, length, 1, 0 if derived else ip, 0 if derived else 40,
                    0 if derived else 1, 0)
    data[mp+24:mp+176] = descriptor(False) + descriptor(True)
    fp = len(data)
    data += b'H10F' + pack('IQII', 1, 48, 1, 0) + pack('IIQQ', 0, int(fanin), mp, 176)
    data[:88] = b'HIC\0' + pack('I10Q', 10, 88+len(variable), fp, 48, 0, 0, 0, 0, 0, 0, 0)
    path.write_bytes(data)


def check(straw, probe, directory):
    from test_v10 import run
    for score, fanin, name in [(False, False, 'large-derived-count'), (True, False, 'large-derived-score'), (True, True, 'large-score-fanin')]:
        path = directory / (name + '.hic')
        fixture(path, score, fanin)
        rows = run([probe, path, 'raw', 'chrA', 'chrB' if fanin else 'chrA', 1000 if fanin else 2]).splitlines()
        if fanin:
            assert rows == ['0 0 s ' + str(int.from_bytes(pack('f', 70000.), 'little'))], rows
        else:
            value = 's 1065353216' if score else 'c 1'
            assert rows == [f'{i} {i+1} {value}' for i in range(70000)]
            # Overlapping cis windows must not double-count after subdivision.
            rows = run([straw, 'observed', 'NONE', path, 'chrA:2:139998', 'chrA:0:139996', 'BP', 2]).splitlines()
            assert len(rows) == 69999
            actual = {tuple(map(int, row.split())) for row in rows}
            assert actual == {(2*i, 2*i+2, 1) for i in range(1, 69998)} | {(2, 0, 1), (139998, 139996, 1)}

if __name__ == '__main__':
    directory = pathlib.Path(sys.argv[1])
    for score, fanin, name in [(False, False, 'large-derived-count'), (True, False, 'large-derived-score'), (True, True, 'large-score-fanin')]:
        fixture(directory / (name + '.hic'), score, fanin)
