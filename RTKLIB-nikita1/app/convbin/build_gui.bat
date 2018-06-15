set prj_dir=build\Win32\Default
"C:\Program Files\CMake\bin\cmake.exe" -G "Visual Studio 15 2017" . -B%prj_dir%
cd %prj_dir%
convbin.sln
