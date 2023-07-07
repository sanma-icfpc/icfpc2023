#!/bin/bash
# This script needs to be run on bash.

echo "building glog"
pushd $(dirname $0)
# do not build in glog/build since it is not listed in glog/.gitignore and thus leads to a dirty state.
mkdir -p glog_build
pushd glog_build
cmake -DBUILD_TESTING:BOOL=OFF -DWITH_GFLAGS:BOOL=OFF ../glog/ -G "Unix Makefiles" || exit 1
make glog || exit 1
popd
popd
