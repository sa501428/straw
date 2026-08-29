# libstraw for .NET

This .NET 8 library uses P/Invoke and `SafeHandle` owners over the stable C ABI.
Install `libstraw` where the platform loader can find it, then reference
`Straw/Straw.csproj`.

```csharp
using var file = new HicStraw.HicFile("sample.hic");
var contacts = file.Records("2", "1", 10_000, normalization: "KR");
```

Queries transfer coordinates and values in three bulk copies. Exact V10 counts
are exposed as `ulong` through `HicFile.Raw`.
