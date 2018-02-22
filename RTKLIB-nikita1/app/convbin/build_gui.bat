set prj_dir=build\Win32\Default
cmake . -B%prj_dir%
cd %prj_dir%
convbin.sln
