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
