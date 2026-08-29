# hic-straw for Ruby

This gem is a thin Ruby FFI adapter over `libstraw`. Install the `ffi` gem and
make `libstraw` discoverable by the platform loader, or set `LIBSTRAW_PATH`.

```ruby
file = HicStraw::File.new("sample.hic")
records = file.records("2", "1", resolution: 10_000, normalization: "KR")
file.close
```

Parsing and decompression remain native, and contact columns cross FFI in bulk.
