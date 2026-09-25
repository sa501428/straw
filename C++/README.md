# Straw C++

This is a C++ implementation of the Straw tool for reading Hi-C data and converting it to a simpler slice format.

This is the only flavor that creates Slice (`.slc`) or HBS (`.hbs.gz`) files.
Those formats are intentionally absent from the C ABI and language wrappers.

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
Half of those fine-resolution windows are guided by a reservoir of occupied
bins from the preceding coarser resolution; the other half remain random so
equal coarse aggregates cannot conceal relocated fine contacts.
Normalization vectors are exhaustive once loaded; expected vectors are
exhaustive at coarse resolutions and deterministically sampled at fine
resolutions. `ALL` is excluded. Missing chromosomes/resolutions or a capability
present in only one file count as differences. A matrix or vector unavailable
in both files is skipped and counted in the summary.

Raw regions from both files are streamed concurrently into a sharded sparse
delta. Exact matches are removed immediately; only non-cancelling cells remain
for tolerance checks and `sum(abs(A-B))` reporting. This avoids retaining two
whole sparse matrices during an exhaustive high-resolution comparison.

```sh
# Practical thorough comparison (the default strategy).
build/straw compare old.v9.hic new.v10.hic

# Compare every stored raw cell and every vector value at every resolution.
build/straw compare old.v9.hic new.v10.hic --all

# More fine-resolution coverage and custom numerical tolerances.
build/straw compare a.hic b.hic --samples 12 --window-bins 512 \
  --abs-tol 1e-5 --rel-tol 1e-5 --seed 42

# Probe an additional/custom normalization name. Repeating --norm replaces
# the default VC, VC_SQRT, KR, and SCALE list.
build/straw compare a.hic b.hic --norm SCALE --norm VC
```

Other controls are `--exhaustive-at BP` and `--max-errors N`. Exit status is 0
when the files match within tolerance, 1 when differences are found, and 2 for
invalid arguments or an unreadable/invalid input. The summary reports the exact
amount of matrix, cell, and vector coverage, so a sampled result is not confused
with a byte-for-byte or fully exhaustive comparison.

### Randomized region sampling

Instead of the fixed 256-bin windows, regions can be drawn at random across
sizes, diagonal distances, and resolutions. Sizes and distances are drawn
log-uniformly so that small and large values get comparable coverage; distances
are drawn in bins rather than base pairs so a coarse resolution does not collapse
every draw onto the diagonal. Everything is seeded, so a run is reproducible.

```sh
# Vary the region size within the normal per-pair sweep.
build/straw compare a.hic b.hic --min-window-bins 8 --max-window-bins 2048

# Sweep intra-chromosomal regions from near-diagonal to 100 Mb apart.
build/straw compare a.hic b.hic --vary-distance --min-distance 0 --max-distance 100000000

# Draw 500 fully random (chromosome pair, resolution, size, distance) regions
# instead of walking every pair, restricted to three random shared resolutions.
build/straw compare a.hic b.hic --random-regions 500 --sample-resolutions 3 \
  --min-window-bins 4 --max-window-bins 1024 --vary-distance --intra-fraction 0.7
```

| Option | Effect |
| --- | --- |
| `--min-window-bins N` / `--max-window-bins N` | draw region widths log-uniformly in that range (default: fixed `--window-bins`) |
| `--vary-distance` | draw the diagonal offset of intra regions log-uniformly |
| `--min-distance BP` / `--max-distance BP` | bound the diagonal offset (either implies `--vary-distance`) |
| `--resolution BP` | restrict to a shared resolution; repeatable |
| `--sample-resolutions N` | randomly keep N of the shared resolutions |
| `--random-regions N` | draw N random regions globally instead of walking every pair |
| `--intra-fraction P` | share of random regions that are intra-chromosomal (default 0.5) |
| `--skip-vectors` | compare regions only, for a pure region-read benchmark |

### Stratified sweep

`--stratified N` draws, **for every resolution**, N regions in each of four
strata: near-diagonal, mid-range, far-from-diagonal, and inter-chromosomal.
Within a stratum the draws cycle through a shuffled list of eligible
chromosomes, so N draws hit N distinct chromosomes (or chromosome pairs) before
repeating. Band edges are `--near-max BP` (default 1 Mb) and `--mid-max BP`
(default 10 Mb). A chromosome too short to host a band is excluded, and a
stratum with no eligible chromosome is reported as skipped rather than silently
dropped.

This is the mode to use for a head-to-head reader benchmark: the report ends
with a per-(resolution, stratum) table of median region-read latency for each
file, the ratio, and which file won. Reads alternate which file goes first, so
the file read second does not pick up a warm-page-cache advantage.

Derived resolutions — those a V10 file stores no matrix for and aggregates from
a finer source at query time — are covered like any other, and are marked
`[derived]` in the timing tables (`[derived:1st]` / `[derived:2nd]` when only one
file computes them on the fly), since they cost the reader extra work.

### Timing metrics

Every query is timed, and a report is printed unless `--no-timing` is given. It
covers region reads, normalization vectors, and expected vectors for each file,
with count, total, mean, p50, p95, max, and throughput, plus breakdowns of region
reads by resolution, by region size, and by diagonal distance. Because each query
re-opens the file and parses its header, the numbers reflect a cold per-query
cost, which is what a caller issuing independent queries actually pays.

```sh
# Benchmark a V9 reader against a V10 reader over a balanced sweep: 25 regions
# per stratum per resolution, sizes from 8 to 512 bins, each region read three
# times, with the per-query rows kept for offline analysis.
build/straw compare file1.v9.hic file1.v10.hic --stratified 25 \
  --min-window-bins 8 --max-window-bins 512 --repeat 3 --timing-csv timings.csv
```

`--repeat N` reads each region N times on both files, timing every read while
comparing only the first; the first read of a region is cold and later ones are
warm, so the spread across repeats separates I/O from decode cost.
`--timing-csv PATH` writes one row per query with columns `type, file, path,
chrom1, chrom2, resolution, norm, x0, x1, y0, y1, width_bins, distance, records,
seconds, iteration`.

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

Positions are numeric BP bin starts computed as `bin_index * resolution`. They
are emitted as stored, without assigning or recording a pair-level coordinate
origin; preprocessors consume the values as supplied. All real cis and trans
chromosome pairs are included, with each stored cell emitted once; the synthetic
`ALL` overview is excluded. Chromosome pairs remain contiguous for the
preprocessors. Cells with zero retained contacts are omitted. Diagnostics go to
stderr.

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
(88-byte header, exact block indexes, and vector indexes), not the older experimental V9-like
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

The reader enforces the fixed high-resolution pyramid: 2 and 5 bp derive from
1 bp, 20 and 50 bp from 10 bp, 200 and 500 bp from 100 bp, and 2 kb from 1 kb.
The 500 kb level must be
materialized. Experimental V10 files that materialized one of those seven virtual
levels are rejected as nonconforming and must be rebuilt.

V10 region values are accepted as supplied without classifying them as zero-based
or one-based. The current `start:end` API uses `end` as an exclusive numeric
boundary and includes bins overlapping that interval. Reversed chromosome queries return coordinates in the requested
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
responses with `Content-Range`. Only candidate blocks and vector chunks intersecting
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
