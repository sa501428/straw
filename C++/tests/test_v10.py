#!/usr/bin/env python3
"""Independent byte-level V10 fixtures; no use of the C++ writer under test."""
import ctypes
import ctypes.util
import functools
import http.server
import pathlib
import struct
import subprocess
import sys
import tempfile
import threading

Z = ctypes.CDLL(ctypes.util.find_library('zstd') or '/opt/homebrew/lib/libzstd.dylib')
Z.ZSTD_compressBound.argtypes = [ctypes.c_size_t]
Z.ZSTD_compressBound.restype = ctypes.c_size_t
Z.ZSTD_compress.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int]
Z.ZSTD_compress.restype = ctypes.c_size_t
Z.ZSTD_isError.argtypes = [ctypes.c_size_t]
Z.ZSTD_isError.restype = ctypes.c_uint

def pack(fmt, *v): return struct.pack('<' + fmt, *v)
def string(s): return s.encode() + b'\0'
def var(n):
    b = bytearray()
    while n >= 128:
        b.append((n & 127) | 128)
        n >>= 7
    b.append(n)
    return bytes(b)
def compress(b):
    out = ctypes.create_string_buffer(Z.ZSTD_compressBound(len(b)))
    n = Z.ZSTD_compress(out, len(out), b, len(b), 3)
    assert not Z.ZSTD_isError(n)
    return out.raw[:n]

def block(rep=0, mode=2, score=False, values=(1, 1, 5), malformed=None, collision=False):
    positions = [0, 4, 5] if collision else [0, 9, 15]
    bitmap = pack('H', sum(1 << p for p in positions))
    if rep == 0: ps = var(positions[0]) + var(positions[1]-positions[0]) + var(positions[2]-positions[1])
    elif rep == 1 or score: ps = bitmap
    else: ps = b''
    scalar = (lambda x: pack('I', x)) if score else var
    if mode == 0: vs = scalar(values[0])
    elif mode == 1:
        vs = scalar(values[0]) + var(1) + var(2) + scalar(values[2])
    elif rep == 2:
        dense = [0] * 16
        for p, v in zip(positions, values): dense[p] = v
        if malformed == 'absent-score': dense[1] = 1
        vs = b''.join(map(scalar, dense))
    else: vs = b''.join(map(scalar, values))
    if malformed == 'overlong': ps = b'\x80\x00' + ps[1:]
    if malformed == 'duplicate': ps = b'\0\0\x0f'
    if malformed == 'zero-count': vs = b'\0' + vs[1:]
    if malformed == 'overflow-varint': vs = b'\xff'*9 + b'\x02' + vs[1:]
    flags = int(rep == 1 or (rep == 2 and score))
    head = pack('5B3x4IQ2I', 1, rep, mode, int(score), flags, 0, 0, 4, 64 if malformed == 'padded-bounds' else 4, 3, len(ps), len(vs))
    assert len(head) == 40
    return head + ps + vs

def fixture(path, rep=0, mode=2, score=False, values=(1, 1, 5), malformed=None,
            transform=0, rotated=False, frag=False, derived=True, collision=False):
    # All data offsets are computed independently from the normative tables.
    chroms = [('chrA', 80), ('chrB', 70)]
    variable = string('test') + pack('I', 2) + string('duplicate') + string('one') + string('duplicate') + string('two')
    variable += pack('I', 2) + b''.join(string(c) + pack('Q', n) for c, n in chroms)
    res = pack('IBBHI', 10, 0, 1, 0, 0xffffffff)
    if derived: res += pack('IBBHI', 20, 1, 1, 0, 0)
    nr = 2 if derived else 1
    variable += pack('I', nr) + res
    variable += pack('I', 1 if frag else 0)
    if frag:
        variable += pack('IBBHI', 1, 0, 1, 0, 0xffffffff)
        for _ in chroms: variable += pack('I', 7) + b''.join(pack('Q', x) for x in range(5, 66, 10))
    variable += pack('I', 1) + string('VC')
    data = bytearray(88) + variable
    ndesc = nr + int(frag)
    matrix_pos = len(data)
    data += b'H10M' + pack('5I', 1, 0, 0, ndesc, 0) + bytes(76 * ndesc)
    logical = block(rep, mode, score, values, malformed, collision)
    directory = var(0) + var(len(logical))
    payload = pack('I', len(directory)) + directory + logical
    frame = compress(payload)
    if malformed == 'concat-frame': frame += compress(b'extra')
    page = b'H10P' + pack('BBHII', 1, 1, 0, len(payload), 1) + frame
    page_pos = len(data)
    data += page
    blob = var(0) + var(len(page)) + var(len(payload))
    index_pos = len(data)
    idx = b'H10I' + pack('5IQ', 1, 1, 64, 1, 0, len(blob))
    idx += pack('4I2Q', 0, 1, 0, 0, page_pos, 0) + blob
    data += idx
    if score:
        floats = [struct.unpack('<f', pack('I', v))[0] for v in values]
        sum_bits = struct.unpack('<Q', pack('d', sum(floats)))[0]
    else: sum_bits = sum(values) & ((1 << 64)-1)
    def descriptor(unit, ri, bin_size, is_derived):
        return pack('4B3IB3xQQII2IQQ2I', unit, int(is_derived), 1, int(score), ri, bin_size,
                    0 if is_derived else 0xffffffff, int(rotated), sum_bits, 1 if is_derived and collision else 3,
                    0x7fc00000, 0x7fc00000, 4, 1 if is_derived else 2,
                    0 if is_derived else index_pos, 0 if is_derived else len(idx),
                    0 if is_derived else 1, 0 if is_derived else 1)
    desc = descriptor(0, 0, 10, False)
    if derived: desc += descriptor(0, 1, 20, True)
    if frag: desc += descriptor(1, 0, 1, False)
    data[matrix_pos + 24:matrix_pos + 24 + len(desc)] = desc
    locators = []
    for kind in range(3):
        entries = []
        for ri in range(nr):
            count = 8 if ri == 0 else 4
            # Derived norm=4 vs source norm=2 catches incorrect reuse.
            value = (2.0 if ri == 0 else 4.0) if kind == 0 else (10.0 if kind == 1 else 20.0)
            words = [struct.unpack('<I', pack('f', value))[0]] * count
            raw = b''.join(pack('I', v) for v in words)
            if transform == 1: raw = b''.join(raw[i::4] for i in range(4))
            if transform == 2: raw = b''.join(pack('I', v ^ (words[i-1] if i else 0)) for i, v in enumerate(words))
            chunk = b'H10V' + pack('BBHII', 1, transform, 0, count*4, count) + compress(raw)
            pos = len(data)
            data += chunk
            chunks = pack('QIBBHQII', 0, count, transform, 1, 0, pos, len(chunk), count*4)
            if kind == 0:
                entry = pack('IIIB3xIIQII', 72, 0, 0, 0, ri, 10*(ri+1), count, 65536, 1) + chunks
            else:
                entry = pack('I', 84 if kind == 2 else 80)
                if kind == 2: entry += pack('I', 0)
                entry += pack('B3xIIQIIII', 0, ri, 10*(ri+1), count, 65536, 1, 1, 0)
                entry += pack('If', 0, 2.0) + chunks
            entries.append(entry)
        index = [b'NVI0', b'EVI0', b'NEVI'][kind] + pack('III', 1, len(entries), 0) + b''.join(entries)
        locators += [len(data), len(index)]
        data += index
    footer_pos = len(data)
    footer = b'H10F' + pack('IQII', 1, 48, 1, 0) + pack('IIQQ', 0, 0, matrix_pos, 24+76*ndesc)
    data += footer
    data[:88] = b'HIC\0' + pack('I10Q', 10, 88+len(variable), footer_pos, len(footer), *locators, 0)
    if malformed == 'bad-footer': struct.pack_into('<Q', data, 16, len(data)+10)
    if malformed == 'unknown-mode': data[matrix_pos+25] = 7
    if malformed == 'invalid-source': struct.pack_into('<I', data, matrix_pos+24+76+12, 1)
    if malformed == 'short': data = data[:-1]
    path.write_bytes(data)
    return data

def run(command, ok=True):
    p = subprocess.run(list(map(str, command)), text=True, capture_output=True)
    assert (p.returncode == 0) == ok, (command, p.returncode, p.stdout, p.stderr)
    return p.stdout

def main():
    straw, probe = sys.argv[1:3]
    with tempfile.TemporaryDirectory() as d:
        path = pathlib.Path(d) / 'test.hic'
        for rep in range(3):
            for mode in ([0, 1, 2] if rep < 2 else [2]):
                values = (1, 1, 1) if mode == 0 else (1, 1, 5)
                fixture(path, rep, mode, values=values)
                lines = run([probe, path, 'raw', 'chrA', 'chrA', 10]).splitlines()
                assert lines == [f'0 0 c {values[0]}', f'1 2 c {values[1]}', f'3 3 c {values[2]}'], lines
                assert len(run([straw, 'observed', 'NONE', path, 'chrA', 'chrA', 'BP', 20]).splitlines()) == 3
        for rep in range(3):
            for mode in ([0, 1, 2] if rep < 2 else [2]):
                values = (0x80000000,) * 3 if mode == 0 else (0x80000000, 0x80000000, 0x7fc01234)
                fixture(path, rep, mode, score=True, values=values, derived=False)
                lines = run([probe, path, 'raw', 'chrA', 'chrA', 10]).splitlines()
                assert lines == [f'0 0 s {values[0]}', f'1 2 s {values[1]}', f'3 3 s {values[2]}'], lines
        for t in range(3):
            fixture(path, transform=t, rotated=True, frag=True)
            assert run([straw, 'observed', 'VC', path, 'chrA', 'chrA', 'BP', 20]).splitlines() == ['0\t0\t0.0625', '0\t20\t0.0625', '20\t20\t0.3125']
            assert '10\t20\t0.025' in run([straw, 'oe', 'VC', path, 'chrA:10:20', 'chrA:20:30', 'BP', 10])
            assert '20\t10\t1' in run([straw, 'observed', 'NONE', path, 'chrA:20:30', 'chrA:10:20', 'BP', 10])
            assert run([straw, 'observed', 'NONE', path, 'chrA:0:10', 'chrA:10:20', 'BP', 10]) == ''
            assert len(run([probe, path, 'raw', 'chrA', 'chrA', 1, 'FRAG']).splitlines()) == 3
            assert run([probe, path, 'callbacks']).splitlines()[-1] == '3'
            run([probe, path, 'meta'])
            matrix = run([straw, 'observed', 'NONE', path, 'chrA:0:40', 'chrA:0:40', 'MATRIX', 10])
            rows = [list(map(float, row.split())) for row in matrix.splitlines()]
            assert rows[1][2] == rows[2][1] == 1
        fixture(path, malformed='padded-bounds')
        assert len(run([probe, path, 'raw', 'chrA', 'chrA', 10]).splitlines()) == 3
        data = bytearray(fixture(path))
        footer = struct.unpack_from('<Q', data, 16)[0]
        matrix = struct.unpack_from('<Q', data, footer+32)[0]
        index = len(data)
        data += b'H10I' + pack('5IQ', 1, 0, 64, 0, 0, 0)
        struct.pack_into('<QQ', data, matrix+24+20, 0, 0)
        struct.pack_into('<QQII', data, matrix+24+52, index, 32, 0, 0)
        path.write_bytes(data)
        assert run([probe, path, 'raw', 'chrA', 'chrA', 10]) == ''
        fixture(path, collision=True)
        assert run([probe, path, 'raw', 'chrA', 'chrA', 20]).strip() == '0 0 c 7'
        float_bits = lambda x: struct.unpack('<I', pack('f', x))[0]
        fixture(path, collision=True, score=True, values=tuple(map(float_bits, (1.25, 0.5, 2.75))))
        assert run([probe, path, 'raw', 'chrA', 'chrA', 20]).strip() == f'0 0 s {float_bits(4.5)}'
        data = bytearray(fixture(path))
        struct.pack_into('<QQ', data, 32, 0, 0)
        path.write_bytes(data)
        assert '0\t0\t10' in run([straw, 'expected', 'VC', path, 'chrA', 'chrA', 'BP', 10])
        run([straw, 'observed', 'VC', path, 'chrA', 'chrA', 'BP', 10], ok=False)
        fixture(path, values=((1 << 53)+1, 1, 5))
        assert f'c {(1 << 53)+1}' in run([probe, path, 'raw', 'chrA', 'chrA', 10])
        assert f'0\t0\t{(1 << 53)+1}' in run([straw, 'observed', 'NONE', path, 'chrA', 'chrA', 'BP', 10])
        for bad in ['overlong', 'duplicate', 'zero-count', 'overflow-varint', 'concat-frame', 'bad-footer', 'unknown-mode', 'invalid-source', 'short']:
            fixture(path, malformed=bad)
            run([straw, 'observed', 'NONE', path, 'chrA', 'chrA', 'BP', 10], ok=False)
        fixture(path, rep=2, score=True, values=(1, 2, 3), malformed='absent-score')
        run([probe, path, 'raw', 'chrA', 'chrA', 10], ok=False)
        fixture(path, values=((1 << 64)-1, 1, 1))
        run([probe, path, 'raw', 'chrA', 'chrA', 20], ok=False)
        if "--http" not in sys.argv:
            print("V10 local reader conformance passed (HTTP test is opt-in)")
            return
        data = fixture(path)
        requests = []
        class Handler(http.server.BaseHTTPRequestHandler):
            def log_message(self, *args): pass
            def do_GET(self):
                value = self.headers.get('Range')
                assert value and value.startswith('bytes=')
                start, end = map(int, value[6:].split('-'))
                requests.append((start, end))
                part = data[start:end+1]
                self.send_response(206)
                self.send_header('Content-Range', f'bytes {start}-{end}/{len(data)}')
                self.send_header('Content-Length', str(len(part)))
                self.end_headers()
                self.wfile.write(part)
        server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            url = f'http://127.0.0.1:{server.server_port}/test.hic'
            assert run([straw, 'observed', 'NONE', url, 'chrA', 'chrA', 'BP', 10]) == run([straw, 'observed', 'NONE', path, 'chrA', 'chrA', 'BP', 10])
            assert all(end-start+1 < len(data) for start, end in requests)
        finally:
            server.shutdown(); server.server_close(); thread.join()
    print('V10 reader: encodings, exact counts, derived, vectors, APIs, corruption, HTTP passed')
if __name__ == '__main__': main()
