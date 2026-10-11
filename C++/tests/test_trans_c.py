"""Run C ABI trans expected/OE queries on independent V10 fixtures."""
import pathlib
import struct
import sys
import tempfile

sys.dont_write_bytecode = True
from test_v10 import fixture, pack, run

with tempfile.TemporaryDirectory() as directory:
    path = pathlib.Path(directory) / 'trans.hic'
    for score in (False, True):
        values = (1, 1, 5)
        if score:
            values = tuple(struct.unpack('<I', pack('f', v))[0] for v in values)
        fixture(path, trans=True, score=score, values=values)
        run([sys.argv[1], path])
print('C ABI trans expected/OE: count/score, raw/normalized, derived, window passed')
