@echo off
rem Runs the Demo. The logic lives in run.sh.
call "%~dp0run-bash.cmd" run.sh %*
exit /b %errorlevel%
