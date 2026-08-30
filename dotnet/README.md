# libstraw for .NET

This .NET 8 library uses P/Invoke and `SafeHandle` owners over the stable C ABI.
The NuGet package carries native runtime assets for `linux-x64`, `osx-x64`,
`osx-arm64`, and `win-x64`; the SDK selects the matching `runtimes/<rid>/native`
library. Source builds can instead install `libstraw` on the loader path.

```csharp
using var file = new HicStraw.HicFile("sample.hic");
var contacts = file.Records("2", "1", 10_000, normalization: "KR");
```

Queries transfer coordinates and values in three bulk copies. Exact V10 counts
are exposed as `ulong` through `HicFile.Raw`.

The current managed surface provides file version/genome, chromosomes,
resolutions, sparse records, and exact V10 raw records. Native failures become
`StrawException` with the stable C status code. Dense matrices, vectors, and
prepared queries remain available in the C ABI but are not yet wrapped here.
Slice and HBS creation are C++-only.
