#!/bin/sh

bash libs/build.sh || exit 1
make -j8 -B -C src dirs solver k3solver tssolver test || exit 1

mkdir -p build || exit 1
cp src/solver build/main || exit 1
