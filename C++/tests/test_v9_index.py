#!/usr/bin/env python3
"""V9 direct-index/fallback parity and HTTP request regression tests."""
import argparse
import ctypes
import http.server
import math
from pathlib import Path
import struct
import sys
import subprocess
import tempfile
import threading
import zlib

BACKEND = "cli"


def pack(fmt, *args):
    return struct.pack('<' + fmt, *args)


def string(value):
    return value.encode() + b'\0'


def fixture(path, direct=True):
    data = bytearray(b'HIC\0' + pack('IQ', 9, 0) + string('hg38'))
    patch = len(data)
    data += bytes(16) + pack('II', 0, 2)
    for name, length in [('chrA', 100), ('chrB', 80)]:
        data += string(name) + pack('Q', length)
    data += pack('III', 1, 10, 0)
    data += bytes(131072 - len(data))
    raw = pack('iii4B5h', 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 8)
    block = zlib.compress(raw)
    block_at = 262144
    matrices = []
    for a, b in [(0, 0), (0, 1)]:
        at = len(data)
        matrix = pack('3i', a, b, 1) + string('BP') + pack('i4f4i', 0, 8, 1, 0, 0, 10, 100, 1, 1)
        matrix += pack('iqi', 0, block_at, len(block))
        data += matrix
        matrices.append((f'{a}_{b}', at, len(matrix)))
    data += bytes(block_at - len(data)) + block
    data += bytes(524288 - len(data))
    master = len(data)
    struct.pack_into('<Q', data, 8, master)
    directory = pack('i', len(matrices))
    for key, at, length in matrices:
        directory += string(key) + pack('qi', at, length)
    raw_begin = master + 8 + len(directory)
    n = 20000
    def expected(normalized):
        return ((string('VC') if normalized else b'') + string('BP') + pack('iQ', 10, n)
                + pack(f'{n}f', *([8 if normalized else 4] * n)) + pack('iif', 1, 0, 2))
    footer = directory + pack('i', 1) + expected(False)
    data += pack('Q', len(footer)) + footer + pack('i', 1) + expected(True)
    nvi = len(data)
    index = bytearray(pack('i', 3))
    payload = bytearray()
    # Include a superseded duplicate key; the last value wins.
    entries = [(0, 4, 11), (1, 4, 9), (0, 2, 11)]
    index_size = 4 + sum(len(string('VC')) + 4 + len(string('BP')) + 4 + 8 + 8 for _ in entries)
    vectors_begin = nvi + index_size + 131072
    for chr_id, value, count in entries:
        vector = pack('Q', count) + pack(f'{count}f', *([value] * count))
        index += string('VC') + pack('i', chr_id) + string('BP') + pack('iQQ', 10, vectors_begin + len(payload), len(vector))
        payload += vector
    data += index + bytes(131072) + payload
    if direct:
        struct.pack_into('<QQ', data, patch, nvi, len(index))
    path.write_bytes(data)
    return dict(master=master, raw_begin=raw_begin, nvi=nvi, index_size=len(index), patch=patch, vectors_begin=vectors_begin)


def query(binary, path, mt='observed', norm='VC', trans=False, ok=True):
    second = 'chrB' if trans else 'chrA'
    if BACKEND == 'c':
        lib = ctypes.CDLL(str(binary))
        ptr = ctypes.c_void_p
        lib.straw_file_open.argtypes = [ctypes.c_char_p, ctypes.POINTER(ptr), ctypes.POINTER(ptr)]
        lib.straw_query_records_simple.argtypes = [ptr] + [ctypes.c_char_p]*5 + [ctypes.c_int32, ctypes.POINTER(ptr), ctypes.POINTER(ptr)]
        lib.straw_records_size.argtypes = [ptr]; lib.straw_records_size.restype = ctypes.c_size_t
        for name, typ in [('x', ctypes.c_int64), ('y', ctypes.c_int64), ('values', ctypes.c_float)]:
            fn = getattr(lib, 'straw_records_' + name)
            fn.argtypes = [ptr]; fn.restype = ctypes.POINTER(typ)
        for name in ('straw_file_close', 'straw_records_free', 'straw_error_free'):
            getattr(lib, name).argtypes = [ptr]
        file, records, error = ptr(), ptr(), ptr()
        try:
            status = lib.straw_file_open(str(path).encode(), ctypes.byref(file), ctypes.byref(error))
            if not status:
                status = lib.straw_query_records_simple(file, mt.encode(), norm.encode(), b'chrA', second.encode(),
                                                       b'BP', 10, ctypes.byref(records), ctypes.byref(error))
            if not ok:
                assert status != 0, 'malformed index was accepted'
                return
            assert status == 0, f'C API error {status}'
            assert lib.straw_records_size(records) == 1
            return (lib.straw_records_x(records)[0], lib.straw_records_y(records)[0], lib.straw_records_values(records)[0])
        finally:
            if error: lib.straw_error_free(error)
            if records: lib.straw_records_free(records)
            if file: lib.straw_file_close(file)
    command = [str(binary), mt, norm, str(path), 'chrA', second, 'BP', '10']
    if BACKEND == 'python':
        # Native calls retain the GIL. Run them in a subprocess so the parent
        # loopback server can serve requests while the extension is executing.
        code = ('import hicstraw,sys; a=sys.argv[1:]; '
                'records=hicstraw.straw(a[0],a[1],a[2],a[3],a[4],"BP",10); '
                '[print(r.binX,r.binY,r.counts) for r in records]')
        command = [sys.executable, '-c', code, mt, norm, str(path), 'chrA', second]
    if BACKEND == 'r':
        code = ('a<-commandArgs(TRUE); d<-strawr::straw(a[2],a[3],a[4],a[5],"BP",10,a[1]); '
                'for(i in seq_len(nrow(d))) cat(d$x[i],d$y[i],format(d$counts[i],digits=15),"\\n")')
        command = [str(binary), '-e', code, mt, norm, str(path), 'chrA', second]
    proc = subprocess.run(command, capture_output=True, text=True, timeout=15)
    if not ok:
        assert proc.returncode != 0, (proc.stdout, proc.stderr)
        return
    assert proc.returncode == 0, proc.stderr
    records = [tuple(map(float, row.split())) for row in proc.stdout.splitlines()]
    assert len(records) == 1, (records, proc.stderr)
    return records[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('binary', type=Path)
    parser.add_argument('--http', action='store_true')
    parser.add_argument('--backend', choices=("cli", "python", "c", "r"), default="cli")
    parser.add_argument('--probe', type=Path)
    args = parser.parse_args()
    global BACKEND
    BACKEND = args.backend
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        direct, fallback = root/'direct.hic', root/'fallback.hic'
        layout = fixture(direct)
        fixture(fallback, False)
        for mt in ('observed', 'oe', 'expected'):
            for norm in ('NONE', 'VC'):
                for trans in (False, True):
                    a = query(args.binary, direct, mt, norm, trans)
                    b = query(args.binary, fallback, mt, norm, trans)
                    assert a == b, (mt, norm, trans, a, b)
        for mt, norm, expected in [('observed', 'VC', 2), ('observed', 'NONE', 8),
                                   ('expected', 'VC', 4), ('oe', 'VC', .5), ('oe', 'NONE', 4)]:
            assert math.isclose(query(args.binary, direct, mt, norm)[2], expected, rel_tol=1e-6)
        # A locator with no length still uses the footer fallback.
        partial = root/'partial.hic'
        content = bytearray(direct.read_bytes())
        struct.pack_into('<Q', content, layout['patch'] + 8, 0)
        partial.write_bytes(content)
        assert query(args.binary, partial) == query(args.binary, direct)
        # A present but malformed locator fails; it must not silently scan elsewhere.
        for name, position, length in [('truncated', layout['nvi'], layout['index_size'] - 1),
                                       ('outside', len(content) + 10, layout['index_size'])]:
            bad = root/f'{name}.hic'
            content = bytearray(direct.read_bytes())
            struct.pack_into('<QQ', content, layout['patch'], position, length)
            bad.write_bytes(content)
            query(args.binary, bad, ok=False)
        if args.probe:
            proc = subprocess.run([str(args.probe), str(direct), 'norm-names'], capture_output=True, text=True, check=True)
            assert set(proc.stdout.splitlines()) == {'NONE', 'VC'}
        if args.http:
            requests = []
            payloads = {'/direct': direct.read_bytes(), '/fallback': fallback.read_bytes()}
            class Handler(http.server.BaseHTTPRequestHandler):
                def log_message(self, *args): pass
                def do_GET(self):
                    data = payloads[self.path]
                    first, last = map(int, self.headers['Range'][6:].split('-'))
                    last = min(last, len(data)-1)
                    requests.append((self.path, first, last))
                    body = data[first:last+1]
                    self.send_response(206)
                    self.send_header('Content-Range', f'bytes {first}-{last}/{len(data)}')
                    self.send_header('Content-Length', str(len(body)))
                    self.end_headers()
                    self.wfile.write(body)
            server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), Handler)
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            try:
                base = f'http://127.0.0.1:{server.server_port}'
                for mt, trans in [('observed', False), ('observed', True), ('oe', True), ('expected', True)]:
                    requests.clear()
                    assert query(args.binary, base+'/direct', mt, trans=trans) == query(args.binary, direct, mt, trans=trans)
                    assert ('/direct', layout['nvi'], layout['nvi'] + layout['index_size'] - 1) in requests, requests
                    assert not any(layout['raw_begin'] <= first < layout['nvi'] or
                                   first <= layout['raw_begin'] + 65536 <= last
                                   for _, first, last in requests), requests
                requests.clear()
                assert query(args.binary, base+'/fallback') == query(args.binary, fallback)
                assert any(layout['raw_begin'] <= first < layout['nvi'] or
                           first <= layout['raw_begin'] + 65536 <= last
                           for _, first, last in requests), requests
                assert query(args.binary, base+'/direct', 'oe') == query(args.binary, direct, 'oe')
                if args.probe:
                    requests.clear()
                    proc = subprocess.run([str(args.probe), base+'/direct', 'norm-names'], capture_output=True, text=True, check=True)
                    assert set(proc.stdout.splitlines()) == {'NONE', 'VC'}
                    assert not any(first >= layout['vectors_begin'] for _, first, _ in requests), requests
            finally:
                server.shutdown(); server.server_close(); thread.join()
    print('V9 direct NVI: values, duplicate keys, missing locator fallback, cis expected/OE, corrupt intervals, HTTP ranges passed')


if __name__ == '__main__':
    main()
