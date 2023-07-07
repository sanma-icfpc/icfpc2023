#!/bin/sh

./build.sh || exit 1
# ICFPC2022時点でのメモ
#   HACK: VS2022 と MSYS2 で、フォルダの位置を合わせる。
#   HACK: C:/home/nodchip/icfpc2022/src/test.exe: error while loading shared libraries: ?: cannot open shared object file: No such file or directory
#         が解消できないので、テストしない。 Windows 側のテストで代用する。
pushd vs/test
../../src/test || exit 1
popd
