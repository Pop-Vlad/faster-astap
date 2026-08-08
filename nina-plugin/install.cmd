@echo off
rem Wrapper so install.ps1 can be double-clicked in Explorer, where the
rem execution policy would otherwise refuse to run it. Arguments are passed
rem straight through, but from a terminal it is simpler to call the script:
rem
rem   powershell -ExecutionPolicy Bypass -File install.ps1 -Uninstall

setlocal
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0install.ps1" %*
set EXITCODE=%ERRORLEVEL%

rem Double-clicked means no arguments and a console of our own that is about to
rem close with the output still in it. Anything with arguments came from a shell
rem that will still be there afterwards, so never block that one.
if "%~1"=="" (
    echo %cmdcmdline% | find /i "%~nx0" >nul
    if not errorlevel 1 pause
)

exit /b %EXITCODE%
