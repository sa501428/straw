# Native MATLAB and Octave binding

The MATLAB flavor calls `libstraw` directly through `straw_mex.cpp`; Python is
not required. Parsing, HTTP I/O, decompression, normalization, and contact
orientation remain in the shared C++ reader.

Build `libstraw`, then configure this directory with CMake and MATLAB:

```sh
cmake -S MATLAB -B MATLAB/build -DCMAKE_PREFIX_PATH=/path/to/libstraw/prefix
cmake --build MATLAB/build
```

Put the resulting `straw_mex` binary and `straw.m` on the MATLAB path.

```matlab
result = straw('KR', 'sample.hic', '1', '1', 'BP', 10000);
[x, y, counts] = straw('NONE', 'sample.hic', '2', '1', 'BP', 100000);
metadata = straw_mex('metadata', 'sample.hic');
raw = straw_mex('raw', 'sample-v10.hic', '1', '1', 'BP', 10000);
```

Sparse coordinates are returned as `int64`, values as `single`, and exact V10
counts as `uint64`.

## Octave

The same source is intended to compile through Octave's MEX compatibility
layer:

```sh
mkoctfile --mex straw_mex.cpp -I/path/to/include -L/path/to/lib -lstraw
```

A separate `.oct` implementation should only be introduced if supported Octave
versions demonstrate an incompatibility or material avoidable copy overhead.
