# hic-straw for Ruby

This gem is a thin Ruby FFI adapter over `libstraw`. Release gems include native
libraries for Linux x64, macOS x64/arm64, and Windows x64. Resolution order is
`LIBSTRAW_PATH`, the bundled platform library, then the system loader.

```ruby
file = HicStraw::File.new("sample.hic")
records = file.records("2", "1", resolution: 10_000, normalization: "KR")
file.close
```

Parsing and decompression remain native, and contact columns cross FFI in bulk.
The Ruby surface includes version/genome metadata, chromosomes, resolutions,
sparse queries, and exact V10 raw counts/scores. Dense/vector/prepared APIs are
not yet wrapped. Slice and HBS creation are C++-only.
