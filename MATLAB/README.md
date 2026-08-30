# Native MATLAB and Octave binding

The MATLAB flavor calls `libstraw` directly through `straw_mex.cpp`; Python is
not required. Parsing, HTTP I/O, decompression, normalization, and contact
orientation remain in the shared C++ reader.

Build `libstraw`, then configure this directory with CMake and MATLAB:

```sh
cmake -S MATLAB -B MATLAB/build -DCMAKE_PREFIX_PATH=/path/to/libstraw/prefix
cmake --build MATLAB/build
```

Release archives contain complete native prefixes under `native/<rid>`. Pass
the matching directory as `-DLIBSTRAW_ROOT=.../native/<rid>` when building the
MEX file; CMake installs both `straw_mex` and `straw.m` together.

Put the resulting `straw_mex` binary and `straw.m` on the MATLAB path.

```matlab
result = straw('KR', 'sample.hic', '1', '1', 'BP', 10000);
[x, y, counts] = straw('NONE', 'sample.hic', '2', '1', 'BP', 100000);
metadata = straw_mex('metadata', 'sample.hic');
raw = straw_mex('raw', 'sample-v10.hic', '1', '1', 'BP', 10000);
```

Sparse coordinates are returned as `int64`, values as `single`, and exact V10
counts as `uint64`.

The binding supports metadata, sparse normalized records, and exact V10 raw
records. Slice and HBS creation are C++-only.

## Octave

The same source is intended to compile through Octave's MEX compatibility
layer:

```sh
mkoctfile --mex straw_mex.cpp -I/path/to/include -L/path/to/lib -lstraw
```

A separate `.oct` implementation should only be introduced if supported Octave
versions demonstrate an incompatibility or material avoidable copy overhead.
