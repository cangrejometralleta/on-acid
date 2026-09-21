@echo off
rem Builds the Package. The logic lives in build.sh.
call "%~dp0run-bash.cmd" build.sh %*
exit /b %errorlevel%
