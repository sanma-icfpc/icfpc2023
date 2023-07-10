IF "%VCINSTALLDIR%" == "" CALL "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" || EXIT /B 1
@ECHO ON
CALL libs\build.bat || EXIT /B 1
REM CALL libs\build_local.bat || EXIT /B 1
msbuild vs\ICFPC2023.sln -t:restore -p:RestorePackagesConfig=true -p:configuration=Release -p:platform=x64 || EXIT /B 1
devenv.com vs\ICFPC2023.sln /rebuild "Release|x64" || EXIT /B 1
