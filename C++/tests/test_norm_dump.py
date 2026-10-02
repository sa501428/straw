"""Exact indexed vector exports, independent of matrices and expected arrays."""
import pathlib
import shlex
import struct
import sys
import tempfile
sys.dont_write_bytecode = True
from test_v10 import fixture, pack, string, compress, run

SPECIAL = [0x3f800001, 0x80000000, 0x7fc01234, 0x7f801234,
           0xffc05678, 0x7f800000, 0xff800000, 0]
WEIRD = 'A #/ "norm"'


def legacy(path, direct=True, big=False):
    data = bytearray(b'HIC\0'+pack('IQ', 9, 0)+string('test'))
    nvi_patch = len(data); data += bytes(16)
    data += pack('II', 0, 3)
    chroms = [('All', 1), ('chrA', 80), ('chrB', 73)]
    for name, n in chroms: data += string(name)+pack('Q', n)
    data += pack('6I', 2, 10, 20, 2, 1, 2)
    for sites in (0, 7, 5): data += pack('I', sites)+b''.join(pack('I', 5+10*i) for i in range(sites))
    # No matrices. Normalized expected deliberately exists to test the fallback
    # NVI locator, but is not exported or needed for stored-vector access.
    struct.pack_into('<Q', data, 8, len(data))
    footer = pack('II', 0, 0)
    data += pack('Q', len(footer))+footer
    data += pack('I', 1)+string('RU')+string('BP')+pack('IQ', 10, 2)
    data += pack('2fIIf', 1, 2, 1, 1, 3)
    source = [('RU', 1, 'BP', 10, [0x40000000]*9),
              ('RU', 0, 'BP', 10, [0x3f800000]),
              (WEIRD, 1, 'BP', 10, SPECIAL),
              ('RU', 1, 'BP', 10, SPECIAL+[0xdeadbeef])]
    wanted = {(WEIRD, 'chrA', 'BP', 10): SPECIAL}
    for norm in ('RU', 'NDSCALE'):
        for ch in (1, 2):
            for unit, bins, size in [('BP', (10, 20), chroms[ch][1]), ('FRAG', (1, 2), (8, 6)[ch-1])]:
                for bin in bins:
                    if norm == 'RU' and ch == 1 and unit == 'BP' and bin == 10:
                        words = SPECIAL+[0xdeadbeef]
                    else:
                        count = size//bin+1  # legacy conventional geometry
                        if norm == 'NDSCALE' and ch == 2 and unit == 'BP' and bin == 20: count = 2
                        if norm == 'NDSCALE' and ch == 1 and unit == 'FRAG' and bin == 2: count = 0
                        words = [0x40000000+i for i in range(count)]
                        source.append((norm, ch, unit, bin, words))
                    wanted[norm, chroms[ch][0], unit, bin] = words
    if big:
        words = (SPECIAL*9000)[:70001]
        source.append(('LONG', 1, 'BP', 10, words))
        wanted['LONG', 'chrA', 'BP', 10] = words
    nvi = len(data)
    index = bytearray(pack('I', len(source))); patches = []
    for norm, ch, unit, bin, words in source:
        index += string(norm)+pack('I', ch)+string(unit)+pack('I', bin)
        patches.append(len(index)); index += pack('QQ', 0, 8+4*len(words))
    payloads = bytearray()
    for patch, (_, _, _, _, words) in zip(patches, source):
        struct.pack_into('<Q', index, patch, nvi+len(index)+len(payloads))
        payloads += pack('Q', len(words))+b''.join(pack('I', v) for v in words)
    data += index+payloads
    if direct: struct.pack_into('<QQ', data, nvi_patch, nvi, len(index))
    path.write_bytes(data)
    return wanted


def read_text(path):
    result = {}; lengths = {}; key = None
    lines = path.read_text().splitlines()
    assert lines[0] == 'HIC_NORM_VECTORS 1'
    for line in lines[1:]:
        fields = shlex.split(line, comments=True)
        if not fields: continue
        if fields[0] == 'vector':
            assert key is None and len(fields) == 5
            key = tuple(fields[1:4])+(int(fields[4]),)
            assert key not in result
            result[key] = []
        elif fields[0] == 'source-length':
            assert key is not None and key not in lengths and not result[key]
            lengths[key] = int(fields[1])
        elif fields[0] == 'end': key = None
        else:
            assert key is not None
            result[key] += [int(field[5:], 16) for field in fields if field.startswith('bits:')]
    assert key is None
    for key, n in lengths.items(): assert len(result[key]) == n
    return result, lengths


def main():
    straw, probe = sys.argv[1:3]
    with tempfile.TemporaryDirectory() as temp:
        root = pathlib.Path(temp)
        run([straw, 'dump-norms', '--help'])
        source = root/'v9.hic'
        for direct in (True, False):
            wanted = legacy(source, direct)
            before = source.read_bytes()
            out = root/f'out-{direct}'
            run([straw, 'dump-norms', source, '--output-dir', out])
            paths = list(out.glob('*.norm.txt')); assert len(paths) == 3
            actual = {}
            for path in paths:
                rows, lengths = read_text(path)
                assert len({k[0] for k in rows}) == 1
                actual.update(rows)
            assert actual == wanted, (actual, wanted)
            assert source.read_bytes() == before
            existing = {p: p.read_bytes() for p in paths}
            run([straw, 'dump-norms', source, '--output-dir', out], ok=False)
            assert existing == {p: p.read_bytes() for p in paths}
            assert not list(out.glob('*.tmp.*'))
            one = root/f'RU-{direct}.txt'
            run([straw, 'dump-norms', source, '--norm', 'RU', '--output', one])
            assert read_text(one)[0] == {k: v for k, v in wanted.items() if k[0] == 'RU'}
            names = run([probe, source, 'norm-names']).splitlines()
            assert set(names) == {'NONE', 'RU', 'NDSCALE', WEIRD}
        for t in range(3):
            path = root/f'v10-{t}.hic'; data = bytearray(fixture(path, transform=t))
            nvi = struct.unpack_from('<Q', data, 32)[0]
            # Replace the first NVI chunk with exact unusual words, independent
            # of the reader. All three transforms must preserve every word.
            raw = b''.join(pack('I', word) for word in SPECIAL)
            if t == 1: raw = b''.join(raw[i::4] for i in range(4))
            if t == 2: raw = b''.join(pack('I', word^(SPECIAL[i-1] if i else 0)) for i, word in enumerate(SPECIAL))
            chunk = b'H10V'+pack('BBHII', 1, t, 0, 32, 8)+compress(raw)
            descriptor = nvi+16+40
            struct.pack_into('<QII', data, descriptor+16, len(data), len(chunk), 32)
            data += chunk; path.write_bytes(data)
            dest = root/f'VC-{t}.txt'
            run([straw, 'dump-norms', path, '--norm', 'VC', '--output', dest])
            rows, lengths = read_text(dest)
            assert rows == {('VC', 'chrA', 'BP', 10): SPECIAL, ('VC', 'chrA', 'BP', 20): [0x40800000]*4}
            assert not lengths
            # A corrupt compressed payload fails without exposing any output.
            data[-len(chunk)+16] = 0; path.write_bytes(data)
            fail = root/f'bad-{t}'
            run([straw, 'dump-norms', path, '--output-dir', fail], ok=False)
            assert not list(fail.iterdir())
        wanted = legacy(source, big=True)
        dest = root/'long.txt'
        run([straw, 'dump-norms', source, '--norm', 'LONG', '--output', dest])
        assert read_text(dest)[0] == {k: v for k, v in wanted.items() if k[0] == 'LONG'}
        legacy(source)
        for args in [('--norm', 'NONE', '--output-dir', root/'bad'),
                     ('--norm', 'absent', '--output-dir', root/'bad'),
                     ('--output', root/'bad.txt'),
                     ('--norm', 'RU', '--norm', 'RU', '--output', root/'bad.txt'),
                     ('--norm', 'RU', '--output', root/'bad.txt', '--output-dir', root/'bad'),
                     ('--output-dir', root/'bad', '--unknown', 'x')]:
            run([straw, 'dump-norms', source, *args], ok=False)
        assert not (root/'bad').exists()
        # Truncated index/vector must be rejected rather than silently outputting
        # zero words or accepting a descriptor count that overreads its extent.
        source.write_bytes(source.read_bytes()[:-1])
        run([straw, 'dump-norms', source, '--output-dir', root/'truncated'], ok=False)
        assert not (root/'truncated').exists()
        if '--http' in sys.argv:
            import http.server
            import threading
            expected = legacy(source, direct=False)
            other = root/'http-v10.hic'; fixture(other, transform=2)
            payloads = {'/v9': source.read_bytes(), '/v10': other.read_bytes()}
            class Handler(http.server.BaseHTTPRequestHandler):
                def log_message(self, *args): pass
                def do_GET(self):
                    data = payloads.get(self.path, payloads['/v9'])
                    first, last = map(int, self.headers['Range'][6:].split('-'))
                    part = data[first:last+1]
                    self.send_response(206)
                    # Wrong offsets must fail, even when the response has the
                    # requested number of bytes (a silent shifted range).
                    offset = 1 if self.path == '/bad-range' else 0
                    self.send_header('Content-Range', f'bytes {first+offset}-{last+offset}/{len(data)}')
                    self.send_header('Content-Length', str(len(part)))
                    self.end_headers(); self.wfile.write(part)
            server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), Handler)
            thread = threading.Thread(target=server.serve_forever, daemon=True); thread.start()
            try:
                base = f'http://127.0.0.1:{server.server_port}'
                dest = root/'http-v9'
                run([straw, 'dump-norms', base+'/v9', '--output-dir', dest])
                actual = {}
                for path in dest.glob('*.norm.txt'): actual.update(read_text(path)[0])
                assert actual == expected
                dest = root/'http-v10.txt'
                run([straw, 'dump-norms', base+'/v10', '--norm', 'VC', '--output', dest])
                assert read_text(dest)[0] == {('VC', 'chrA', 'BP', 10): [0x40000000]*8,
                                              ('VC', 'chrA', 'BP', 20): [0x40800000]*4}
                run([straw, 'dump-norms', base+'/bad-range', '--output-dir', root/'bad-range'], ok=False)
                assert not (root/'bad-range').exists()
            finally:
                server.shutdown(); server.server_close(); thread.join()
    print('Normalization dump tests passed')

if __name__ == '__main__': main()
