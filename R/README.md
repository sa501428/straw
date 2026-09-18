# strawr

`strawr` reads legacy V6-V9 files and consolidated-format V10 files. V10
support includes derived resolutions, Zstandard-compressed blocks and vectors,
BP and FRAG queries, normalization, and expected-value matrices. Numeric V10
region values are consumed as supplied without requiring or inferring a
coordinate origin; the current `start:end` API treats `end` as an exclusive
numeric boundary.
Straw is a fast implementation of reading/dump for .hic files

The R flavor exposes `.hic` metadata and sparse contact queries. It does not
create Slice or HBS files; those outputs are C++-only.

## Installation
```R
remotes::install_github("aidenlab/straw/R")
```

## Usage
```R
# <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>
hic.data.frame <- strawr::straw("KR", "/path/to/file.hic", "11", "11", "BP", 10000)
```
