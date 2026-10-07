@echo off
setlocal
pushd "%~dp0"
python -B -X utf8 -m lagcheck wizard
set "lagcheck_exit=%errorlevel%"
popd
exit /b %lagcheck_exit%
