# libstraw C API

This directory contains the stable C-compatible ABI for the C++ straw reader.
It is an adapter over the maintained C++ implementation, not an independent C
parser.

All public declarations are in `include/straw/straw.h`. Native result objects
own their arrays; accessor pointers remain valid until the matching
`straw_*_free` function is called. Errors are explicit objects and C++
exceptions never cross the ABI.

See `ABI.md` for the ownership, orientation, threading, and compatibility
contract.

Build from `C++/`:

```sh
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

The resulting shared library is `libstraw`. See `examples/query.c` for a small
pure-C consumer.

## Functionality

`straw_file_open` retains file metadata and exposes chromosomes, resolutions,
normalizations, attributes, vectors, record counts, sparse records, dense
matrices, and exact V10 raw `uint64_t` counts. `straw_query_prepare` now owns a
parsed native reader plus prepared matrix/zoom state; repeated `window` and
`regions` calls do not reconstruct the reader or footer state, and the query
remains valid if its original `straw_file_t` is closed.

Failures are classified from typed native exceptions, never by matching error
message text. Callers receive a stable `straw_status_t` and optional owned
`straw_error_t`.

Tagged releases provide install prefixes containing `include/`, `lib/` (and
`bin/` on Windows), CMake package metadata, and `pkg-config` metadata. Set
`CMAKE_PREFIX_PATH` to an unpacked prefix or use `pkg-config --cflags --libs
libstraw`. Slice and HBS are intentionally not exposed by this ABI.
