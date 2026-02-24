# Straw C++

Fast C++ implementation for reading `.hic` files and extracting contact matrix data. Supports region queries, genome-wide binary slice dumps, and plain-text dumps of raw contact data.

---

## Dependencies

| Dependency | Purpose |
|---|---|
| CMake >= 3.13 | Build system |
| libcurl | Remote `.hic` file access over HTTP |
| zlib | Block decompression |
| zstd | Block decompression (v10+ files) |
| C++14 compiler | Language standard |

### Install dependencies

**macOS (Homebrew):**
```bash
brew install cmake curl zlib zstd
```

**Ubuntu/Debian:**
```bash
sudo apt install cmake libcurl4-openssl-dev zlib1g-dev libzstd-dev
```

**RHEL/Fedora:**
```bash
sudo dnf install cmake libcurl-devel zlib-devel libzstd-devel
```

---

## Build

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

The binary is placed at `build/straw`.

---

## Usage

### 1. Region query (sparse text output)

Extract contact records for a chromosome pair or sub-region:

```
straw [observed|oe|expected] <NONE|VC|VC_SQRT|KR> <hicFile> \
    <chr1>[:x1:x2] <chr2>[:y1:y2] <BP|FRAG|MATRIX> <binsize>
```

Output format (tab-separated):
```
<binX>  <binY>  <value>
```

**Examples:**
```bash
# Whole chromosome, no normalization
straw observed NONE input.hic chr1 chr1 BP 10000

# Sub-region with KR normalization
straw observed KR input.hic chr1:0:5000000 chr2:0:5000000 BP 5000

# Output as dense matrix
straw observed NONE input.hic chr1 chr1 MATRIX 1000000
```

---

### 2. `dump` — genome-wide binary slice file

Dump all (or filtered) chromosome pairs into a binary `.slc` file:

```
straw dump <observed|oe|expected> <NONE|VC|VC_SQRT|KR> <hicFile> \
    <BP|FRAG> <binsize> <outputFile> <compressed> [filter]
```

| Argument | Values | Description |
|---|---|---|
| `compressed` | `1` / `true` / `compressed` | gzip-compress the output |
| `filter` | _(omit)_ | All contacts |
| | `-intra` | Intra-chromosomal only |
| | `-intra-short` | Intra within 5 Mb |
| | `-intra-long` | Intra beyond 5 Mb |
| | `-inter` | Inter-chromosomal only |

**Examples:**
```bash
# All contacts, compressed, 10 kb resolution
straw dump observed NONE input.hic BP 10000 output.slc 1

# Intra-chromosomal only, uncompressed
straw dump observed NONE input.hic BP 10000 output.slc 0 -intra

# Short-range contacts only
straw dump observed NONE input.hic BP 5000 output.slc 1 -intra-short
```

#### Slice file format (`.slc`)

| Field | Type | Description |
|---|---|---|
| Magic | 8 bytes | `"HICSLICE"` |
| Resolution | `int32` | Bin size in bp |
| numChromosomes | `int32` | Number of chromosomes |
| Chromosome table | repeated | `int32` name length, name bytes, `int16` key |
| Records | repeated | `int16` chr1Key, `int32` binX, `int16` chr2Key, `int32` binY, `float` value |

The compressed variant wraps the entire file in gzip.

---

### 3. `dump-text` — genome-wide plain-text dump

Dump all chromosome pairs at a given resolution as raw (NONE normalization) plain text to stdout:

```
straw dump-text <hicFile> <resolution>
```

Output format (space-separated):
```
<chr1> <pos1> <chr2> <pos2> <value>
```

Zero, NaN, and Inf values are skipped. Only the upper triangle is emitted (`chr1.index <= chr2.index`).

**Examples:**
```bash
# Dump entire genome at 1 Mb resolution
straw dump-text input.hic 1000000

# Save to file
straw dump-text input.hic 10000 > contacts.txt

# Pipe into another tool
straw dump-text input.hic 25000 | awk '$5 > 10'
```

---

## Notes

- For region queries, chromosome indices are automatically swapped so that `chr1.index <= chr2.index` as required by the `.hic` format. You do not need to specify chromosomes in any particular order.
- The `dump-text` and `dump` commands both read directly from the compressed block data without loading the full matrix into memory, making them suitable for large files.
- HTTP URLs are supported wherever a file path is accepted (requires libcurl).
