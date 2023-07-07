#!/bin/sh

#./build.sh || exit 1
# HACK: VS2022 と WSL/MSYS2 で、フォルダの位置を合わせる。
pushd vs/test
rm -f test_result_linux.xml
../../src/test --gtest_output=xml:test_result_linux.xml || exit 1
popd
