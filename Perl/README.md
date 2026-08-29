# Hic::Straw for Perl

`Hic::Straw` uses `FFI::Platypus` to call the stable `libstraw` C ABI. Install
the native library and set `LIBSTRAW_PATH` if it is not on the loader path.

```perl
my $file = Hic::Straw->open('sample.hic');
my $records = $file->records('2', '1', resolution => 10_000);
$file->close;
```

The Perl layer only converts bulk native arrays; it does not parse `.hic` data.
