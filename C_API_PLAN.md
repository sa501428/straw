# Cross-Language C API Plan

## Objective

Expose the maintained C++ straw reader through a small, stable C-compatible
ABI, then build thin language-specific wrappers around that ABI. Parsing,
indexed I/O, HTTP range access, decompression, normalization, and expected-value
handling will remain in the C++ implementation.

This is a C API over the C++ implementation, not a separate C rewrite of the
`.hic` parser.

The official cross-language API will continue to exclude the subsample and
comparison modes. Those features may remain available through the C++ command
line program, but they will not constrain or expand the stable library ABI.

## Scope

The first stable API should cover the reusable reader functionality that is
currently exposed through C++ entry points:

- Local `.hic` files and HTTP(S) URLs
- V6 through V10 files
- File format version and genome identifier
- Chromosome metadata
- BP and FRAG resolutions
- Available normalization methods and file attributes
- Sparse region queries
- Dense matrix queries
- `observed`, `oe`, and `expected` values
- Arbitrary normalization names, including `NONE`
- Normalization vectors
- Expected-value vectors
- Batched region queries
- Genome-wide and per-chromosome record counts
- Exact V10 raw integer counts and stored float scores
- Consistent contact orientation in the order requested by the user for every
  supported `.hic` version and every maintained straw flavor
- Existing slice and HBS export operations, after the read-only API is stable

The following are explicitly out of scope:

- Subsampling
- File comparison
- CLI argument parsing
- CLI text formatting
- Independent parsers implemented in each target language

## Cross-Cutting Contact Orientation Bug Fix

Before freezing the C ABI, define and enforce one contact-orientation contract
across all supported `.hic` versions and straw flavors:

- The first returned coordinate (`x` or `binX`) always corresponds to the first
  chromosome or region supplied by the user.
- The second returned coordinate (`y` or `binY`) always corresponds to the
  second chromosome or region supplied by the user.
- Internal canonical chromosome ordering is an implementation detail and must
  never change the orientation visible to the caller.
- Reversing an interchromosomal query from `(A, B)` to `(B, A)` returns the same
  contacts and values with the coordinates transposed.
- Dense results for the reversed query are the transpose of the forward result,
  subject to the requested row and column windows.
- For cis queries with asymmetric row and column windows, a stored canonical
  contact is reflected only when necessary to satisfy the requested axes. A
  contact must not be dropped, returned on the wrong axes, or emitted twice.
- Normalization vectors, expected-value calculations, region filters, and
  bounds checks must be applied to the logical user axes, not accidentally to
  the internally canonicalized axes.

The fix must cover:

- V6, V7, V8, V9, and V10 inputs
- Local files and HTTP(S) sources
- BP and FRAG queries
- Sparse, dense, streamed, batched-region, normalized, OE, expected, and exact
  raw query paths where each operation is supported
- The C++ library and CLI
- The pybind11 Python flavor
- `strawr`
- The current MATLAB access path until it is retired, and the new MATLAB/Octave
  native binding
- The new C ABI and all wrappers built on it
- `straw-rust`, including its native legacy reader and V10 bridge

Java and JavaScript straw implementations live in separate repositories. Audit
them against the same contract and either submit corresponding fixes and tests
there or record that they already conform. The shared behavior should not be
called complete while an actively maintained official flavor is known to return
the opposite orientation.

Implement the correction in the lowest shared reader/query layer available.
Language wrappers should not each invent their own flip heuristic. Wrappers may
assert and test the contract, but the native query result they consume should
already be oriented to the user's request.

### Orientation regression fixtures

Add deliberately asymmetric fixtures or generated test files for every format
family. They should contain distinguishable chromosome lengths, chromosome
indices, regions, coordinates, contact values, and normalization-vector values
so that an incorrect swap cannot pass accidentally.

At minimum, test:

- A forward interchromosomal query where file chromosome order matches request
  order
- The same query with chromosome order reversed
- Different row and column subregions, also queried in reverse
- Cis queries with disjoint or partially overlapping asymmetric windows
- Sparse and dense parity, including the dense transpose relationship
- `observed NONE`, normalized observed, OE, and expected results
- Exact V10 raw counts and V10 stored scores
- BP and FRAG units
- Batched regions containing both empty and nonempty results
- Local and HTTP access

For a forward result containing `(x, y, value)`, the equivalent reversed query
must contain `(y, x, value)`, with the same multiplicity and without duplicates.
Tests should compare contact multisets when output ordering is not part of the
documented API.

## Repository Layout

Add a top-level `C/` directory for the public ABI and its adapter. The C++
parser and query implementation remains under `C++/`.

```text
straw/
|-- C++/
|   |-- straw.cpp
|   |-- straw.h
|   |-- straw_v10.cpp
|   `-- ...
|-- C/
|   |-- include/
|   |   `-- straw/
|   |       `-- straw.h
|   |-- src/
|   |   `-- straw_c.cpp
|   |-- tests/
|   |   `-- test_straw_c.c
|   |-- examples/
|   |   `-- query.c
|   |-- CMakeLists.txt
|   `-- README.md
|-- MATLAB/
|-- Julia/
|-- Ruby/
|-- Perl/
|-- dotnet/
`-- Go/                 # Added only if the Go wrapper is implemented
```

`C/include/straw/straw.h` must be valid C. `C/src/straw_c.cpp` will implement
that API in C++ and translate between C data structures and the existing C++
reader.

The installed public header should be:

```text
include/straw/straw.h
```

The installed shared library should use the normal platform name for
`libstraw`, such as `libstraw.so`, `libstraw.dylib`, or `straw.dll`.

## Build Architecture

Refactor the build into reusable targets instead of compiling implementation
sources separately into every consumer.

```text
straw_core       Internal C++ parser and query engine
|-- straw        Existing command-line executable
`-- libstraw     Shared/static library exporting the public C ABI
```

The C adapter, CLI, Python extension, and other native consumers must link to a
defined library target. They should not include `straw.cpp` directly or build
private copies of the parser.

The CMake work should provide:

- Shared and optional static library builds
- Correct public and private include directories
- Hidden symbol visibility by default
- An export macro for the public C symbols
- Install rules for the library and header
- CMake package configuration
- `pkg-config` metadata
- Platform library versioning and `SONAME`/install-name handling
- Windows DLL export and import support
- Continued construction of the existing `straw` executable

CLI-only sources such as `subsample.cpp` and `compare.cpp` should not become
part of the public library ABI.

## C ABI Design

### General rules

- No C++ class, exception, `std::string`, `std::vector`, callback object, or
  standard-library allocator may cross the ABI boundary.
- All exported symbols use the `straw_` prefix.
- Public handles are opaque.
- All fallible operations return a stable status code.
- C++ exceptions are caught at the ABI boundary.
- Error categories are machine-readable; wrappers must not parse error text.
- Strings crossing the API are UTF-8.
- All output arguments are initialized to a safe value on failure.
- Every native allocation has one documented owner and one matching release
  function.
- Public coordinates use fixed-width types and do not depend on C `long` size.

### Opaque handles

The initial API should use handles similar to:

```c
typedef struct straw_file straw_file_t;
typedef struct straw_query straw_query_t;
typedef struct straw_error straw_error_t;
typedef struct straw_records straw_records_t;
typedef struct straw_raw_records straw_raw_records_t;
typedef struct straw_dense_matrix straw_dense_matrix_t;
typedef struct straw_vector straw_vector_t;
typedef struct straw_metadata straw_metadata_t;
```

An open file handle retains parsed file metadata. A prepared query handle
retains the chromosome pair, matrix type, unit, resolution, and normalization
needed for repeated region queries.

### Error handling

Define stable error categories, including:

- Success
- Invalid argument
- I/O or HTTP failure
- Unsupported file version
- Corrupt or truncated file
- Chromosome not found
- Resolution or unit not available
- Normalization or expected vector not available
- Unsupported operation
- Allocation failure
- Internal error

A representative operation is:

```c
straw_status_t straw_file_open(
    const char *path_or_url,
    straw_file_t **out_file,
    straw_error_t **out_error);

void straw_file_close(straw_file_t *file);
void straw_error_free(straw_error_t *error);
```

The error object should expose its status, a diagnostic message, and, where
available, contextual information such as the input path or unavailable
capability. Passing `NULL` as the error output should be allowed when the caller
only needs the status code.

### ABI versioning

Provide:

- A compile-time ABI version macro
- `straw_abi_version()`
- `straw_library_version()`
- A documented ABI compatibility policy
- Structure-size fields for public option structures that may grow
- Reserved fields where appropriate
- An exported-symbol allowlist checked in CI

Compatible additions retain the ABI major version. Breaking changes require a
new ABI major version rather than silently changing an existing structure or
symbol.

## Memory and Result Model

### Sparse records

Return sparse results in bulk, using a structure-of-arrays representation:

```text
x:      int64[count]
y:      int64[count]
values: float32[count]
```

One result object owns all three arrays. Accessors return borrowed pointers that
remain valid until `straw_records_free()` is called. Callers must not free the
individual pointers.

The public API should use 64-bit coordinates even though the existing legacy
`contactRecord` uses 32-bit coordinates. The stable ABI should not preserve that
known limitation.

The normal sparse result represents the established float-valued straw query
surface. It is suitable for observed, normalized, expected, and OE values.

### Exact V10 raw records

Keep the exact raw API separate from the general float-valued query API. V10
raw records may contain either exact `uint64_t` counts or stored float scores.
The raw result must preserve that distinction and must not route exact counts
through `float` or `double`.

The exact API should expose:

- `uint64_t` bin indices
- Whether the returned values are counts or scores
- Exact `uint64_t` counts
- Stored `float32` scores

Legacy files store raw values with legacy float precision. The API must document
that this source precision cannot be recovered as exact integers.

### Dense matrices

Return a single contiguous `float32` allocation with explicit:

- Row count
- Column count
- Storage order
- Leading dimension, if needed
- Query origin or interval
- Resolution

Choose and document one canonical storage order before ABI v1 is frozen.
Benchmark MATLAB conversion before making the decision because MATLAB is
column-major while most C-facing consumers expect row-major storage.

### Vectors, strings, and metadata

- Normalization and expected vectors use contiguous `double` arrays.
- Resolutions use contiguous `int32_t` arrays.
- Chromosome lengths use `int64_t`.
- String collections are owned by one result object and accessed by index.
- All result objects have a type-specific `straw_*_free()` function.
- Borrowed metadata pointers remain valid for the documented lifetime of their
  file or result handle.

### Batched regions

Expose batched queries as one concatenated result:

```text
region_offsets: size_t[region_count + 1]
x:              int64[total_records]
y:              int64[total_records]
values:         float32[total_records]
```

This prevents a foreign-function call per contact or per small region and is a
high-value path for MATLAB, Julia, and .NET workloads.

## C API Surface

### Library information

- ABI version
- Library version
- Optional build and feature information

### File and metadata operations

- Open and close a local path or URL
- File format version
- Genome identifier
- Chromosome list
- BP resolutions
- FRAG resolutions
- Available normalization names
- File attributes
- Normalization availability by chromosome, unit, and resolution

### Query operations

- One-shot sparse query
- Prepared query creation and destruction
- Sparse window query through a prepared query
- Dense matrix query
- Batched region query
- Normalization vector retrieval
- Expected vector retrieval
- Exact V10 raw query
- Genome-wide record count
- Per-chromosome record counts

Use enums for closed sets such as matrix type, unit, and storage order.
Normalization remains a UTF-8 string because normalization names are extensible.

### Export operations

After the read-only API is stable, expose existing non-subsampling export
operations through a versioned options structure:

- Slice output
- HBS output
- Existing intra/inter filters
- Existing compression and output options

Export should retain transactional output behavior where it is already
available. Subsampling and comparison remain outside this surface.

## Implementation Milestones

### Milestone 0: Contact orientation correctness

This correctness milestone precedes the reusable core and ABI freeze so the new
API does not standardize inconsistent legacy behavior.

1. Write the user-requested-axis contract as a C++ regression test.
2. Audit every V6-V9 sparse, streaming, dense, multi-region, vector, and count
   path for internal chromosome canonicalization.
3. Audit the equivalent V10 general and exact-raw paths.
4. Fix region swapping, returned-coordinate swapping, dense transposition, and
   normalization-axis selection in the shared reader layer.
5. Port the regression cases to Python, R, MATLAB, and `straw-rust`.
6. Audit the separately maintained Java and JavaScript flavors and track any
   required fixes in their repositories.
7. Document the orientation contract in user-facing API documentation.

Acceptance criteria:

- Every supported format version returns `x` for the first user query axis and
  `y` for the second.
- Reversed sparse queries contain exactly the transposed contacts and equal
  values, without loss or duplication.
- Reversed dense queries have the documented transpose relationship.
- Asymmetric cis windows return the correct canonical or reflected contact once.
- Normalized and OE queries use normalization data for the correct logical axis.
- All maintained flavors pass equivalent orientation regression cases.
- No language wrapper contains a version-specific workaround that disagrees
  with the shared contract.

### Milestone 1: Reusable C++ core

1. Create `straw_core` from the maintained reader sources.
2. Link the existing CLI against `straw_core`.
3. Separate CLI-only code from reusable reader code.
4. Stop directly including `.cpp` implementation files in native bindings.
5. Preserve existing CLI behavior and test results.

Acceptance criteria:

- The existing C++ executable builds and passes its current tests.
- At least one small C++ consumer links to `straw_core` without compiling its
  own copy of the reader.
- Subsample and comparison code is absent from the public library target.

### Milestone 2: C directory and minimal ABI

1. Add the proposed `C/` directory.
2. Add a pure-C public header.
3. Implement library version, file open/close, error handling, chromosome
   metadata, and resolution metadata.
4. Add a C example and a smoke test compiled as C11.
5. Add install rules for the shared library and header.

Acceptance criteria:

- A C compiler can compile the public header without C++ extensions.
- A C program can open a fixture, inspect chromosomes and resolutions, handle
  an invalid input, and release every allocation.
- No C++ symbol or type is required by the consumer.

### Milestone 3: Core query ABI

1. Implement sparse bulk queries.
2. Implement normalization and expected vectors.
3. Implement dense matrices.
4. Implement record-count operations.
5. Implement prepared query handles.
6. Implement batched region queries.
7. Implement exact V10 raw records.

Acceptance criteria:

- C results match the C++ API on V6 through V10 fixtures.
- Exact V10 counts above `2^53` remain exact.
- Repeated prepared queries do not reopen and reparse the file unnecessarily.
- Empty results, missing capabilities, and corrupt inputs return stable statuses.

### Milestone 4: Hardening and packaging

1. Add sanitizer testing for leaks and invalid memory access.
2. Add independent-handle concurrency tests.
3. Check the exported symbol allowlist.
4. Add Linux, macOS, and Windows builds.
5. Produce installable native artifacts.
6. Publish the ABI and ownership documentation.

Acceptance criteria:

- Repeated open/query/free and error paths are clean under ASan and LSan.
- Independent handles work concurrently as documented.
- Installed examples build using only installed headers and library metadata.
- Native artifacts can be consumed without a source-tree-relative build.

### Milestone 5: Export API

1. Add slice and HBS export options to the C ABI.
2. Preserve existing filtering and output semantics.
3. Test failures without leaving a partially published output file where
   transactional behavior is promised.

This milestone can follow the first language bindings if it would otherwise
delay the high-value read/query API.

## Shared Conformance Suite

Create language-neutral fixtures and expected results. The pure-C test driver is
the reference test of the ABI, and each wrapper runs the same behavioral cases.

Coverage should include:

- Local V6 through V10 files
- HTTP byte-range access
- Materialized and derived V10 resolutions
- BP and FRAG units
- Sparse, dense, observed, normalized, OE, and expected queries
- Normalization and expected vectors
- Reversed chromosome order
- User-axis contact orientation for V6 through V10 and every matrix/value mode
- Asymmetric interchromosomal and cis windows, including reflected cis contacts
- Region boundary and empty-result behavior
- Batched regions
- Missing chromosome, resolution, normalization, and expected data
- Corrupt or truncated files
- Exact integer counts above `2^53`
- Stored V10 float scores
- All allocation and free paths
- Repeated open and close cycles
- Independent-handle concurrency

Wrapper results should be checked for array length, shape, coordinates, values,
error category, and exact integers where exactness is promised.

## Language Bindings

### Julia: high priority

Use Julia's native `ccall`/`@ccall` support directly.

Provide:

- An `HicFile` object with explicit `close` and a finalizer
- Metadata functions
- Sparse and dense queries
- Normalization and expected vectors
- Batched queries
- Exact raw values using `UInt64`
- Artifact-managed native libraries where feasible

The safe default should copy results into Julia-owned memory. An advanced
zero-copy view may be offered only while the native result owner is kept alive.

### MATLAB

Replace the Python-mediated path with a direct native MEX gateway linked to
`libstraw`.

Provide:

- A function compatible with the familiar straw query form
- MATLAB-native sparse and dense outputs
- Metadata functions
- `uint64` exact raw counts
- Optional open-file and prepared-query handles for repeated queries
- Deterministic cleanup using `mexAtExit`
- Appropriate `mexLock`/`mexUnlock` handling for retained handles

The wrapper should never require a Python installation.

### Octave

First compile and test the MATLAB MEX gateway through Octave's MEX compatibility
layer. Keep one gateway source with small compatibility macros where possible.

Add a native `.oct` wrapper only if tests or benchmarks show that:

- Required MEX APIs are unavailable,
- Handle cleanup cannot be made reliable,
- Data conversion introduces material avoidable copying, or
- Supported Octave platforms cannot load the MEX build consistently.

MATLAB and Octave should run the same conformance cases.

### C#/.NET

Use P/Invoke over the bulk C API.

Provide:

- `SafeHandle` implementations for native owners
- `IDisposable` public objects
- UTF-8 string handling
- Checked native-length conversions
- Bulk copying or safe spans over native memory
- Runtime-specific NuGet assets for supported systems

Do not perform one P/Invoke call per record.

### Ruby

Begin with Ruby FFI over the C API. Copy bulk arrays into Ruby-owned objects at
the public boundary and provide deterministic or block-scoped cleanup.

Add a small native extension only if profiling shows that FFI conversion or Ruby
object construction is a material bottleneck. Parsing and decompression remain
in C++ in either case.

### Perl

Begin with `FFI::Platypus`, with destructors for native handles and bulk array
conversion. Use XS only if packaging or measured conversion overhead warrants
it. Do not implement a pure-Perl `.hic` parser.

### Go

Add a cgo wrapper only after there is demonstrated user demand and a clear
maintenance owner. Copy native bulk results into Go slices and do not retain C
pointers after the owning result is freed.

Avoid callbacks in the initial Go surface.

## Existing Binding Consolidation

After the ABI has remained stable for at least one release:

1. Change the Python extension to link to library targets instead of including
   implementation `.cpp` files.
2. Consider using the C ABI for Python's simple functions while retaining
   pybind11 where it provides useful object ergonomics.
3. Migrate `strawr` away from its duplicated legacy parser toward `libstraw`.
4. Evaluate changing the Rust V10 bridge to consume the stable library ABI
   instead of compiling C++ implementation sources directly.

These migrations should not block the initial C ABI or first-wave bindings, but
they are necessary to reach the long-term goal of one maintained parser.

## Recommended Delivery Order

1. Define the user-requested-axis contract and add failing orientation fixtures.
2. Fix contact flipping across V6 through V10 and all existing straw flavors.
3. Create the reusable `straw_core` CMake target.
4. Add `C/` and a minimal pure-C metadata API.
5. Freeze the ABI v1 status, error, versioning, ownership, and orientation
   conventions.
6. Implement sparse records and vector access.
7. Implement dense, prepared, batched, exact-raw, and record-count operations.
8. Complete the pure-C conformance suite and native packaging.
9. Implement Julia and MATLAB in parallel as the first binding wave.
10. Validate the MATLAB MEX gateway under Octave.
11. Implement .NET, followed by Ruby and Perl.
12. Add slice/HBS export if it was deferred from ABI v1.
13. Consolidate the existing Python and R native implementations.
14. Implement Go only in response to demonstrated demand.

## Definition of Done

The cross-language effort is complete when:

- `libstraw` exposes a documented and versioned C ABI.
- The public header compiles as C without C++ dependencies.
- No C++ type or exception crosses the ABI.
- Every allocation has an unambiguous owner and release function.
- Large results cross language boundaries in bulk.
- V6 through V10 results match the maintained C++ implementation.
- Every maintained straw flavor returns contacts in the user-requested axis
  order for V6 through V10, including reversed chromosome requests and
  asymmetric cis windows.
- Exact V10 raw integer counts remain exact.
- MATLAB no longer requires Python.
- The MEX gateway works in Octave, or a measured reason for a `.oct` wrapper is
  documented.
- Julia, MATLAB/Octave, .NET, Ruby, and Perl pass the shared conformance suite.
- Subsample and comparison remain outside the official cross-language API.
- Supported platforms have documented native installation or packaging paths.
