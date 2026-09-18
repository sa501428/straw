# libstraw ABI and ownership contract

`libstraw` ABI version 1 is declared by `straw/straw.h`. The header is C11 and
may also be included from C++. All exported functions use the `straw_` prefix;
no C++ type, exception, or allocator crosses the boundary.

## Ownership

- `straw_file_open` creates a file handle owned by the caller. Release it with
  `straw_file_close`.
- `straw_query_prepare` creates a prepared query owned by the caller. Release
  it with `straw_query_close`. It retains its own parsed reader and matrix/zoom
  state, remains valid after the originating file handle is closed, and reuses
  that state for every window or batch call.
- Query, vector, dense, raw, and batch functions create one opaque result owner.
  Release it with the matching type-specific `straw_*_free` function.
- Array accessors return borrowed pointers. They remain valid until their result
  owner is freed and must never be freed individually.
- File metadata strings are borrowed from the file handle and remain valid
  until `straw_file_close`.
- A non-null error returned through `out_error` is owned by the caller and must
  be released with `straw_error_free`.
- Every output owner is set to null before a fallible operation begins.
- Every function accepts a null `out_error` when the caller only needs the
  status code.

## Result representation

Sparse and batched contacts use structure-of-arrays storage. Coordinates are
signed 64-bit genomic or fragment coordinates and values are float32, matching
the established general straw query API. Dense matrices are contiguous
row-major float32 arrays.

The separate V10 raw API returns uint64 bin indices and parallel count, score,
and kind arrays. Count records retain their exact uint64 values; score records
retain the stored float32 value. Exact raw access is rejected for legacy files,
whose source values already use legacy floating-point storage.

Numeric query coordinates are consumed as supplied. The ABI does not require,
infer, convert, or record whether a caller's values originated in a zero-based
or one-based convention. Stored/raw bin indices remain ordinary array indices.

## Contact orientation

The first coordinate always corresponds to the first chromosome or region in
the user query, and the second coordinate corresponds to the second. Reversing
an interchromosomal query transposes its coordinates. Reflected cis contacts are
returned once in the requested row/column orientation. Language wrappers must
not add their own version-specific flip rules.

## Threading

Independent file, query, and result handles may be used concurrently. A result
must not be freed while another thread is reading its borrowed pointers. ABI v1
does not promise that simultaneous operations on the same prepared query handle
share mutable cached state.

## Errors

Native failures cross the ABI as `straw_status_t` plus an optional owned error
message. Status selection is based on typed C++ exceptions and explicit
capability checks; message text is descriptive only and is never parsed to
determine the status.

## Compatibility

The high 16 bits of `straw_abi_version()` are the ABI major version and the low
16 bits are the minor version. Compatible additions increment the minor
version. Existing symbols, enum values, ownership rules, and structure prefixes
remain unchanged within ABI major version 1. A breaking change requires a new
ABI major version.
