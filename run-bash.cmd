@echo off
rem Finds a bash and hands it the script named by the first argument.
rem This is the one helper every shim calls; no shim searches on its own.
rem
rem The Order of the Search, widest toolbox first:
rem   1. Git for Windows, which Shares the Windows PATH where Python lives.
rem   2. A bash on the PATH that is not the WSL launcher in System32.
rem   3. WSL, with the path Translated and said out loud.
rem   4. Nothing: an address to install from, and a non-zero code.
setlocal enabledelayedexpansion

set "SCRIPT=%~1"
if "%SCRIPT%"=="" (
    echo [X] wanted the name of a script to run. 1>&2
    echo     run-bash.cmd build.sh 1>&2
    exit /b 2
)
shift

set "ARGS="
:collect
if "%~1"=="" goto ready
set "ARGS=!ARGS! %1"
shift
goto collect

:ready
cd /d "%~dp0"

rem 1. Git for Windows.
for %%G in ("%ProgramFiles%\Git\bin\bash.exe" "%ProgramFiles(x86)%\Git\bin\bash.exe" "%LocalAppData%\Programs\Git\bin\bash.exe") do (
    if exist "%%~G" (
        "%%~G" "%SCRIPT%"!ARGS!
        exit /b !errorlevel!
    )
)

rem 2. A bash on the PATH, the System32 launcher Discarded.
for /f "delims=" %%B in ('where bash.exe 2^>nul') do (
    echo %%B | find /i "\System32\" >nul
    if errorlevel 1 (
        "%%B" "%SCRIPT%"!ARGS!
        exit /b !errorlevel!
    )
)

rem 3. WSL, path Translated and said out loud.
where wsl.exe >nul 2>&1
if not errorlevel 1 (
    echo -- no Windows bash found; running %SCRIPT% under WSL.
    for /f "delims=" %%P in ('wsl.exe wslpath -a "%CD%" 2^>nul') do (
        wsl.exe --cd "%%P" -- ./%SCRIPT%!ARGS!
        exit /b !errorlevel!
    )
)

rem 4. Nothing.
echo [X] wanted a bash to run %SCRIPT% with, and found none. 1>&2
echo     install Git for Windows: https://git-scm.com/download/win 1>&2
exit /b 1
