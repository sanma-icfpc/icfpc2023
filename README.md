# icfpc2023
ICFPC2023 sanma team repository

## Prerequisites

* Linux (Ubuntu 22.04 LTS (not 18.04 LTS which lacks g++-10)) / WSL2 on Windows 11
  * `sudo apt install build-essential git-core cmake g++-10 clang-15 libtool autoconf texinfo libopencv-dev`
* Windows 10/11
  * Visual Studio 2022 (17.6.3)
  * CMake 3.24.1 (add to PATH)

## Build (Linux / WSL)

```
git clone --recursive https://github.com/sanma-icfpc/icfpc2023.git
# if you forgot to clone with --recursive, try: git submodule update --init
cd icfpc2023
bash build.sh
```

## Build (Visual Studio)

```
git clone --recursive https://github.com/sanma-icfpc/icfpc2023.git
cd icfpc2023
libs/build.bat # to build external libraries (Debug and Release)
start vs/ICFPC2023.sln
select Release;x64 or Debug;x64
Build Solution
```

## Run Tests

```
./test # run all tests
./test --gtest_filter=TestExample.* # run specific tests.
```
