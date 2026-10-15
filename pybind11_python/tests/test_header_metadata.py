"""Header inspection must not need the normalization index or vector payloads."""
import pathlib
import struct
import tempfile
import unittest

import hicstraw


class HeaderMetadataTest(unittest.TestCase):
    def test_v9_header_without_normalization_payloads(self):
        # The index locator points beyond this header-only fixture. A metadata
        # constructor that enumerates normalizations fails on this file.
        data = b'HIC\0' + struct.pack('<IQ', 9, 1_000_000) + b'hg38\0'
        data += struct.pack('<QQII', 2_000_000, 100, 0, 2)
        for name, length in ((b'All', 1000), (b'chr1', 1_000_000)):
            data += name + b'\0' + struct.pack('<Q', length)
        data += struct.pack('<4I', 2, 1000, 5000, 0)
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / 'header.hic'
            path.write_bytes(data)
            hic = hicstraw.HiCFile(str(path))
            self.assertEqual(hic.getGenomeID(), 'hg38')
            self.assertEqual(hic.getResolutions(), [1000, 5000])
            self.assertEqual([(c.name, c.length) for c in hic.getChromosomes()],
                             [('All', 1000), ('chr1', 1_000_000)])


if __name__ == '__main__':
    unittest.main()
