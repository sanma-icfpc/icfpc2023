@echo =================================================
@echo building glog
@echo =================================================
IF "%VCINSTALLDIR%" == "" CALL "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" || EXIT /B 1
pushd %~dp0

if NOT EXIST opencv\x64\vc16\bin\opencv_world470.dll (
mkdir opencv
pushd opencv
curl -o opencv-4.7.0-windows.exe -L https://github.com/opencv/opencv/releases/download/4.7.0/opencv-4.7.0-windows.exe
opencv-4.7.0-windows.exe -o"extract" -y
move extract\opencv\build\include .
move extract\opencv\build\x64 .
popd
)

mkdir glog_vsbuild
pushd glog_vsbuild
cmake -DBUILD_TESTING:BOOL=OFF -DWITH_GFLAGS:BOOL=OFF ..\glog -G "Visual Studio 17 2022" || exit /b 1
MSBuild glog.sln -target:glog -p:Configuration=Release;Platform=x64 || exit /b 1
MSBuild glog.sln -target:glog -p:Configuration=Debug;Platform=x64 || exit /b 1
popd
popd
