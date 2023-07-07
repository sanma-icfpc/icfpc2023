CALL build.bat
PUSHD vs\test
..\x64\Release\test.exe || EXIT /B 1
POPD
