# Hic::Straw for Perl

`Hic::Straw` uses `FFI::Platypus` to call the stable `libstraw` C ABI. Release
archives include `lib/Hic/Straw/native/<rid>`; loading prefers
`LIBSTRAW_PATH`, then the
bundled library, then the system loader.

```perl
my $file = Hic::Straw->open('sample.hic');
my $records = $file->records('2', '1', resolution => 10_000);
$file->close;
```

The Perl layer only converts bulk native arrays; it does not parse `.hic` data.
It currently exposes version/genome, chromosomes, resolutions, and sparse
records. Raw/dense/vector/prepared APIs are not yet wrapped. Slice and HBS
creation are C++-only.
