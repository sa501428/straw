"""Run after building the Python extension in place; requires NumPy."""
import pathlib
import struct
import sys
import tempfile

root = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / 'pybind11_python'))
sys.path.insert(0, str(root / 'C++/tests'))
sys.dont_write_bytecode = True
from test_v10 import fixture, pack
import hicstraw

with tempfile.TemporaryDirectory() as directory:
    path = pathlib.Path(directory) / 'trans.hic'
    for score in (False, True):
        values = (1, 1, 5)
        if score:
            values = tuple(struct.unpack('<I', pack('f', v))[0] for v in values)
        fixture(path, trans=True, score=score, values=values, frag=True)
        for resolution, mean, product in ((10, 7/56, 6), (20, 7/16, 24)):
            for norm in ('NONE', 'VC'):
                for mt in ('oe', 'expected'):
                    for a, b in (('chrA', 'chrB'), ('chrB', 'chrA')):
                        records = hicstraw.straw(mt, norm, str(path), a, b, 'BP', resolution)
                        expected = [mean]*3 if mt == 'expected' else [v/mean/(product if norm == 'VC' else 1) for v in (1,1,5)]
                        assert len(records) == 3
                        assert all(abs(r.counts-v) < 1e-5 for r, v in zip(records, expected))
                zoom = hicstraw.HiCFile(str(path)).getMatrixZoomData('chrB', 'chrA', 'oe', norm, 'BP', resolution)
                expected = 1/mean/(product if norm == 'VC' else 1)
                window = zoom.getRecords(0, 10, 0, 10)
                assert len(window) == 1 and abs(window[0].counts-expected) < 1e-5
                assert abs(float(zoom.getRecordsAsMatrix(0, 10, 0, 10)[0,0])-expected) < 1e-5
                assert abs(float(hicstraw.strawAsMatrix('oe', norm, str(path),
                    'chrB:0:10', 'chrA:0:10', 'BP', resolution)[0,0])-expected) < 1e-5
        frag = hicstraw.straw('oe', 'NONE', str(path), 'chrB', 'chrA', 'FRAG', 1)
        assert len(frag) == 3 and abs(frag[0].counts-64/7) < 1e-5
print('Python trans expected/OE: count/score, raw/normalized, derived, reversed, window, dense, FRAG passed')
