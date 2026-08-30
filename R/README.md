# strawr

`strawr` reads legacy V6-V9 files and consolidated-format V10 files. V10
support includes derived resolutions, Zstandard-compressed blocks and vectors,
BP and FRAG queries, normalization, and expected-value matrices. V10 interval
ends are exclusive, as defined by the V10 format specification.
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
