# Straw.jl

Thin Julia bindings for the stable `libstraw` C ABI. Parsing and decompression
remain in C++; results cross the boundary in bulk arrays.

Release archives include platform libraries under `native/<rid>`. Loading tries
`LIBSTRAW_PATH` first, then the bundled library, `deps/usr/lib`, and finally the
system loader. Source users may unpack a matching `libstraw` release and set
`LIBSTRAW_PATH` to the full library path.

```julia
using Straw

file = HicFile("sample.hic")
contacts = records(file, "1", "1"; normalization="KR", resolution=10_000)
matrix = dense(file, "1:0:100000", "1:0:100000"; resolution=10_000)
close(file)
```

The module exposes metadata, sparse and dense records, stored normalization and
expected vectors, exact V10 raw counts/scores, and retained prepared single or
batched queries. Results own Julia arrays after each native call. Slice and HBS
creation are C++-only.
