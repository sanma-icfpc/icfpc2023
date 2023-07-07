CALL build.bat
PUSHD vs\test
if exist test_result_windows.xml del test_result_windows.xml
..\x64\Release\test.exe --gtest_output=xml:test_result_windows.xml || EXIT /B 1
POPD
