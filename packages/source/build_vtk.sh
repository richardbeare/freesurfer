#!/usr/bin/env bash

set -e

[ "$#" != "1" ] && echo "error: usage: build.sh <prefix>" && exit 1
INSTALL_DIR="$1"

cd vtk

cmake . \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_FLAGS="-DGLX_GLXEXT_LEGACY" \
  -DCMAKE_CXX_FLAGS="-DGLX_GLXEXT_LEGACY" \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
  -DVTK_WRAP_TCL=ON

make -j8
make install
