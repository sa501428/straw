# Straw C++

This is a C++ implementation of the Straw tool for reading Hi-C data and converting it to a simpler slice format.

## Description:
The tool provides functionality to read .hic files and extract contact matrices at specific resolutions. It can output the data in a simplified binary slice format that is more efficient for downstream processing at single high-resolutions.

## Installation:
1. Requires CMake 3.13 or higher
2. Requires libcurl, zlib, and zstd development libraries, plus pkg-config
3. Clone the repository
4. Create a build directory: `mkdir build`
5. Enter build directory: `cd build`
6. Run cmake: `cmake ..`
7. Build: `make`

## Usage:
The main executable 'straw' supports three modes:
1. Standard mode:
`straw [observed/oe/expected] <NONE/VC/VC_SQRT/KR> <hicFile> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG/MATRIX> <binsize>`
2. Dump mode (creates slice file):
`straw dump <observed/oe/expected> <NONE/VC/VC_SQRT/KR> <hicFile> <BP/FRAG> <binsize> <outputFile>`
3. Subsample mode (prints weighted short text):
`straw subsample <hicFile> <--fraction P|--contacts N> [--resolution BP] [--seed N] [--output output.hbs.gz]`
4. Compare mode (validates two `.hic` files):
`straw compare <first.hic> <second.hic> [options]`

## Examples:
1. Extract specific region:
`straw observed NONE input.hic chr1:0:1000000 chr2:0:1000000 BP 10000`
2. Create slice file at 10kb resolution:
`straw dump observed NONE input.hic BP 10000 output.slc`

## Compare two Hi-C files

`compare` checks chromosome names and lengths, advertised BP resolutions, raw
nonzero contacts, normalization vectors, raw expected vectors, and normalized
expected vectors. It understands both V6–V9 and V10, and uses the exact V10 raw
API so integer counts are not rounded to float before comparison.

By default every chromosome pair is checked exhaustively at shared resolutions
of 100 kb and coarser (and always at the coarsest shared resolution). At finer
resolutions it checks four deterministic 256-by-256-bin windows for every pair.
Vectors are exhaustive at coarse resolutions and deterministically sampled at
fine resolutions. `ALL` is excluded. Missing chromosomes/resolutions or a
capability present in only one file count as differences.

```sh
# Practical thorough comparison (the default strategy).
build/straw compare old.v9.hic new.v10.hic

# Compare every stored raw cell and every vector value at every resolution.
build/straw compare old.v9.hic new.v10.hic --all

# More fine-resolution coverage and custom numerical tolerances.
build/straw compare a.hic b.hic --samples 12 --window-bins 512 \
  --abs-tol 1e-5 --rel-tol 1e-5 --seed 42

# Probe an additional/custom normalization name. Repeating --norm replaces
# the default VC, VC_SQRT, and KR list.
build/straw compare a.hic b.hic --norm SCALE --norm VC
```

Other controls are `--exhaustive-at BP` and `--max-errors N`. Exit status is 0
when the files match within tolerance, 1 when differences are found, and 2 for
invalid arguments or an unreadable/invalid input. The summary reports the exact
amount of matrix, cell, and vector coverage, so a sampled result is not confused
with a byte-for-byte or fully exhaustive comparison.

## Subsample to weighted short text

Export raw `observed NONE` contacts from V6–V9 or V10 to stdout:

```sh
# Retain each contact with probability 0.1.
build/straw subsample input.hic --fraction 0.1 --seed 42 > sampled.short

# Calculate p = 100000000 / total raw contacts at the coarsest BP resolution.
build/straw subsample input.hic --contacts 100000000 --seed 42 > sampled.short

# Export all counts, optionally choosing a particular available BP resolution.
build/straw subsample input.hic --fraction 1 --resolution 1000 > full.short
```

Output has no header and uses five tab-separated columns:

```text
chr1    pos1    chr2    pos2    count
```

Positions are zero-based BP bin starts. All real cis and trans chromosome pairs
are included, with each stored cell emitted once; the synthetic `ALL` overview
is excluded. Chromosome pairs remain contiguous for the preprocessors. Cells
with zero retained contacts are omitted. Diagnostics go to stderr.

The default is the finest BP resolution with nonzero real-chromosome counts.
Empty resolutions are skipped automatically, including ALL-only resolutions
advertised by some converted files. `--resolution` uses the requested resolution
without this fallback. `--contacts` always sums counts at the largest advertised
BP bin size, across all real chromosome pairs, to set its probability. This is
an **expected target**, not an exact retained total; source totals can also
differ between resolutions because of filtering or legacy float rounding.
A target above that coarse total is rejected. Zero is allowed for either mode.

For each cell with raw count `n`, the retained count is distributed as
`Binomial(n, p)`: each underlying contact is assessed independently. A cell
containing two contacts can retain zero, one, or both. Large counts use a
recursive beta/order-statistic binomial sampler, avoiding a loop over reads or
a normal/Poisson approximation. Seed defaults to `1`. The same input, options,
and executable give reproducible output; C++ standard-library differences can
change the sequence between platforms/builds.

Counts must be finite, nonnegative integers. Fractional/negative/nonfinite raw
scores fail instead of being rounded or silently discarded. V10 integer counts
are read and printed as `uint64_t`; legacy float counts retain only the precision
already present in the source. Like other stdout tools, an error can leave
partial output: check the exit status before using the text to rebuild a file.

### Rebuild with hictools-c

For an export at 1000 BP, for example:

```sh
hic_pre -f short -r 1000,5000,10000 sampled.short sampled.v9.hic chrom.sizes
hic_v10 pre -f short -r 1000,5000,10000 sampled.short sampled.v10.hic chrom.sizes
```

Use chromosome names and lengths matching the source. Original read positions,
strands, and restriction fragments cannot be recovered from a binned matrix.
Rebuild at the export resolution or coarser multiples of it, and recompute any
desired normalizations. Currently, hictools-c parses each text weight as float32,
so individual weights above `2^24` may round when rebuilt, even though this
exporter prints V10 integer counts exactly. V9 matrix storage also uses float32.

Subsampling tests run under CTest. To additionally exercise short-format round
trips through both hictools-c builders:

```sh
python3 tests/test_subsample.py build/straw /path/to/hic_pre /path/to/hic_v10
```

## Compressed binary short output (.hbs.gz)

Both `subsample` and `dump` can write HBS, a compact binary alternative to
five-column text:

```sh
build/straw subsample input.hic --fraction 0.1 --seed 42 --output sampled.hbs.gz
build/straw subsample input.hic --contacts 100000000 -o sampled.hbs.gz
build/straw dump observed NONE input.hic BP 1000 full.hbs.gz
build/straw dump observed NONE input.hic BP 1000 trans.hbs.gz -inter
```

`--output` / `-o` selects a named `.hbs.gz` file; subsampling without it continues
to print text to stdout. `dump` selects HBS by the output suffix. Its legacy
compression argument may be omitted for HBS, or supplied for compatibility;
HBS is always gzip-compressed. Existing dump filters are supported. Other output
suffixes retain the existing slice behavior. HBS only accepts `observed NONE`
and `BP`, since it stores integer raw counts.

The header stores resolution, chromosome names, and lengths. Records use uint16
chromosome IDs, uint32 bin indices, and uint16 counts: **14 bytes** before gzip.
Counts of 65,535 or more use an escape followed by uint64, making those records
22 bytes. There is no integer precision loss. See [HBS_FORMAT.md](https://github.com/sa501428/hic-format)
for the versioned byte layout.

Updated hictools-c builders read HBS directly, with auto-detection or `-f hbs`:

```sh
hic_pre -r 1000,5000,10000 sampled.hbs.gz sampled.v9.hic chrom.sizes
hic_v10 pre -r 1000,5000,10000 sampled.hbs.gz sampled.v10.hic chrom.sizes
```

The builders match chromosomes by name, validate lengths, and require output
resolutions to be multiples of the export resolution. V10 preserves exact uint64
counts through this input path; V9 retains its float32 limit. Both currently
require BP positions to fit signed int32. HBS output is staged and published only
after success, so an export error leaves an existing destination untouched.

## Slice Format:
The slice format (.slc) is a binary format that contains:
1. Magic string "HICSLICE"
2. Resolution (int32)
3. Number of chromosomes (int32)
4. Chromosome mapping (name lengths, names, and keys)
5. Contact records (chr1Key, binX, chr2Key, binY, value)

## Reading Slice Files:
A C++ reader is provided in the slice_reader directory. It provides methods to:
1. Read basic file information (resolution, chromosomes)
2. Read all contact records
3. Read records for specific chromosome pairs
The reader automatically handles the chromosome ordering convention.

## Notes:
The simplified slice format is only intended for repeated analysis on a high resolution slice of the matrix. Otherwise, the original hic file format is more efficient.

## Bug Reports or Feature Requests:
For bug reports or feature requests, please open an issue on the repository.
## V10 files

V10 support follows the consolidated `hic-format/HiCFormatV10.md` specification
(88-byte header and page/vector indexes), not the older experimental V9-like
block extension. V6–V9 continue to use their existing read path. The new format
implementation is isolated in `straw_v10.cpp`, `straw_v10.h`, and `v10_binary.h`;
CMake links it automatically. Custom builds must compile `straw_v10.cpp` alongside
`straw.cpp` and link zstd, zlib, curl, and threads.

The existing commands automatically detect V10:

```sh
straw observed NONE input.v10.hic chr1:0:100000 chr1:0:100000 BP 1000
straw observed VC input.v10.hic chr1 chr1 BP 1000
straw oe VC input.v10.hic chr1 chr1 BP 1000
```

Supported V10 features include materialized and derived resolutions, mandatory
rotated cis grids, rectangular trans grids, sparse/bitmap/dense blocks, all value modes, integer
counts and float scores, BP/FRAG metadata, and all three compressed vector
transforms. Derived matrices aggregate raw source cells before applying the
**target** resolution's normalization and expected vectors. Missing capabilities
and corrupt records raise errors; a missing chromosome pair is an empty matrix.
V10 expected/OE queries are defined for cis matrices only.

The reader enforces the fixed high-resolution pyramid: 20 and 50 bp derive from
10 bp, 200 and 500 bp from 100 bp, and 2 kb from 1 kb. The 500 kb level must be
materialized. Experimental V10 files that materialized one of those five virtual
levels are rejected as nonconforming and must be rebuilt.

V10 region ends are **exclusive**. Queries include bins overlapping the requested
interval. Reversed chromosome queries return coordinates in the requested
chromosome order. Cis sparse queries emit each canonical cell once, transposing
it when only its reflected position intersects the window; `MATRIX` output fills
both symmetric entries. FRAG locations use fragment coordinates, and reported
coordinates are fragment-bin starts. These conventions do not change the legacy
V6–V9 path.

### Exact counts and the C++ API

The raw `observed NONE` CLI prints `uint64_t` counts exactly, including values
above 2^53. The existing `contactRecord` APIs and slice files still use `float`,
so their returned counts can round. Use the separate exact API when integer
precision or score bit preservation matters:

```cpp
#include "straw_v10.h"

straw_v10::File file("input.v10.hic");
file.streamRaw("chr1:0:100000", "chr1:0:100000", "BP", 1000,
    [](const straw_v10::Record& record) {
        // binX/binY are bin indices, not base-pair positions.
        // isScore=false: count is the exact uint64_t value.
        // isScore=true: score retains its stored float32 bits.
    });
```

`File::raw` accepts half-open bin ranges directly. The existing metadata,
streaming, region, normalization, dense-matrix, record-count, and slice entry
points also dispatch to V10. Legacy APIs reject coordinate overflow rather than
wrapping a 32-bit coordinate. Parsing, decompression, full vectors, and dense
output have explicit allocation limits (512 MiB per record/allocation class).
Large whole-matrix queries can still require substantial memory, particularly
when aggregating a derived resolution.

### Remote reads and tests

HTTP(S) reading uses byte ranges and requires a server returning exact `206`
responses with `Content-Range`. Only candidate pages and vector chunks intersecting
the query are fetched; a server ignoring Range is rejected instead of downloading
the whole file. Local format tests use independently generated binary fixtures:

```sh
cmake -S . -B build
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

Python tests load the system zstd library through `ctypes`. The HTTP integration
test starts a loopback server and is opt-in:

```sh
cmake -S . -B build -DSTRAW_TEST_HTTP=ON
ctest --test-dir build -R v10_http_ranges --output-on-failure
```
