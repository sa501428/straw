# Straw: rapidly stream data from .hic files
Straw is a family of readers for indexed access to `.hic` contact data. This
repository contains the C++, stable C ABI, Python, R, Julia, .NET, Ruby, Perl,
and MATLAB/Octave flavors. Java and Rust live in sibling repositories.

- Java: https://github.com/aidenlab/java-straw
- Rust: https://github.com/aidenlab/straw-rust

The maintained Python flavor wraps C++ with pybind11. The former pure-Python
implementation is archived at https://github.com/aidenlab/pystraw.

## Flavor status

| Flavor | Implementation and primary surface |
|---|---|
| C++ | Native CLI/library; sparse/dense queries, metadata, vectors, exact V10 raw data, compare, subsample, Slice and HBS |
| C (`libstraw`) | Stable shared-library ABI; owned results, typed status/errors, metadata, sparse/dense/raw/vector queries, retained prepared queries |
| Python / R | In-process C++ bindings for query and metadata workflows |
| Julia / .NET / Ruby / Perl / MATLAB | Thin bindings over release-built `libstraw`; see each directory README for its exposed subset |
| Java | Pure Java V6–V10 reader; caller-oriented axes and `BigInteger` exact V10 counts |
| Rust | Native legacy reader plus maintained V10 bridge; queries, metadata, vectors, exact raw records, retained queries |

Slice (`.slc`) and HBS creation are deliberately C++-only. They are not part of
the C ABI or any other language flavor. JavaScript is outside the scope of this
release work.

## Native releases

Tags matching `v*` run `.github/workflows/native-release.yml`. It tests and
publishes `libstraw` prefixes for Linux x64, macOS x64/arm64, and Windows x64,
then creates native-enabled Julia, .NET, Ruby, Perl, and MATLAB package
archives. A manual workflow run builds the same downloadable CI artifacts
without creating a GitHub release.

A Jupyter notebook example of using straw can be found here: https://aidenlab.gitbook.io/juicebox/accessing-raw-data

## Install straw for python

Use `pip install hic-straw`. 
If you want to build from the source code, you must have pybind11 installed. 
Clone the library and `cd` into the `straw/` directory. Then `pip install ./pybind11_python`.

## Compile straw for C++

```bash
cmake -S C++ -B C++/build -DCMAKE_BUILD_TYPE=Release
cmake --build C++/build --parallel
```

Development headers for cURL, zlib, and zstd are required.
Please see [the wiki](https://github.com/aidenlab/straw/wiki) for more documentation.

For questions, please use
[the Google Group](https://groups.google.com/forum/#!forum/3d-genomics).

Ongoing development work is carried out by <a href="http://mshamim.com">Muhammad S. Shamim</a>.
Past contributors include <a href="http://www.cherniavsky.net/neva/">Neva C. Durand</a> and many others.

If you use this tool in your work, please cite 

**Neva C. Durand, James T. Robinson, Muhammad S. Shamim, Ido Machol, Jill P. Mesirov, Eric S. Lander, and Erez Lieberman Aiden. "Juicebox provides a visualization system for Hi-C contact maps with unlimited zoom." Cell Systems 3(1), 2016.**
