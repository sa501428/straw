#!/usr/bin/env bash
set -euo pipefail

assets=${1:?usage: stage_native_bindings.sh NATIVE_ASSETS_DIR}
root=$(cd "$(dirname "$0")/../.." && pwd)

stage_one() {
  artifact=$1
  rid=$2
  library=$3
  source="$assets/$artifact/stage"
  test -d "$source"

  if test -f "$source/bin/$library"; then
    native="$source/bin/$library"
  else
    case $library in
      *.dylib) pattern='libstraw*.dylib' ;;
      *.so) pattern='libstraw.so*' ;;
      *) pattern="$library" ;;
    esac
    native=$(find "$source/lib" -maxdepth 1 -type f -name "$pattern" | head -n 1)
  fi
  test -n "$native"

  for destination in \
    "$root/Julia/native/$rid/$library" \
    "$root/dotnet/Straw/runtimes/$rid/native/$library" \
    "$root/Ruby/lib/hic_straw/native/$rid/$library" \
    "$root/Perl/lib/Hic/Straw/native/$rid/$library"; do
    mkdir -p "$(dirname "$destination")"
    cp "$native" "$destination"
  done
  mkdir -p "$root/MATLAB/native/$rid"
  cp -R "$source/." "$root/MATLAB/native/$rid/"
}

stage_one libstraw-linux-x64 linux-x64 libstraw.so
stage_one libstraw-osx-x64 osx-x64 libstraw.dylib
stage_one libstraw-osx-arm64 osx-arm64 libstraw.dylib
stage_one libstraw-win-x64 win-x64 straw.dll
