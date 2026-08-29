# Straw.jl

Thin Julia bindings for the stable `libstraw` C ABI. Parsing and decompression
remain in C++; results cross the boundary in bulk arrays.

Install `libstraw` and make it discoverable by the system loader, or set
`LIBSTRAW_PATH` to its absolute path before loading the module.

```julia
using Straw

file = HicFile("sample.hic")
contacts = records(file, "1", "1"; normalization="KR", resolution=10_000)
matrix = dense(file, "1:0:100000", "1:0:100000"; resolution=10_000)
close(file)
```
