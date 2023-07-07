@echo =================================================
@echo building glog
@echo =================================================
IF "%VCINSTALLDIR%" == "" CALL "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" || EXIT /B 1
pushd %~dp0
mkdir glog_vsbuild
pushd glog_vsbuild
cmake -DBUILD_TESTING:BOOL=OFF -DWITH_GFLAGS:BOOL=OFF ..\glog -G "Visual Studio 17 2022" || exit /b 1
MSBuild glog.sln -target:glog -p:Configuration=Release;Platform=x64 || exit /b 1
MSBuild glog.sln -target:glog -p:Configuration=Debug;Platform=x64 || exit /b 1
popd
popd
