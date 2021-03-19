#!/bin/bash
usage="$(basename "$0") [-h] [-r] -- Build script for DronePoseLib

Parameters:
    -h, --help     Show this help message
    -r, --release  Build with optimization flags enabled (for timing)"

set -euo pipefail

release=
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -r|--release) release=1 ;;
        -h|--help) echo "$usage"; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

build_dir=build
build_type=Debug
if [ -n "$release" ]; then
  build_dir=build_release
  build_type=Release
fi

# Display information
echo "Configuring build for DronePoseLib"
echo "  - Build type: $build_type"
echo "  - Build directory: $build_dir"

# Change directory depending on whether it is debug or release mode
mkdir -p $build_dir && cd $build_dir

# Configure
cmake -DCMAKE_BUILD_TYPE=$build_type ..

# Build (for Make on Unix equivalent to `make -j $(nproc)`)
if [[ "$OSTYPE" == "darwin"* ]]; then
  # Mac OSX does not have nproc
  cmake --build . --config $build_type -- -j $(sysctl -n hw.logicalcpu)
else
  cmake --build . --config $build_type -- -j $(nproc)
fi

# Run test script or timing depending on build type
if [ -n "$release" ]; then
  ./example
else
  ctest -j $(nproc) --output-on-failure
fi
