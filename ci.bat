CALL build.bat || EXIT /B 1
PUSHD vs
if exist test\test_result_windows.xml del test\test_result_windows.xml
x64\Release\test.exe --gtest_output=xml:test\test_result_windows.xml || EXIT /B 1
POPD
