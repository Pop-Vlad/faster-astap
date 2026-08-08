@echo off
rem The other side of install.cmd: removes the plugin folder from
rem %localappdata%\NINA\Plugins\3.0.0\. Arguments are passed straight through,
rem so double-clicking removes the plugin and leaves the index cache: several
rem gigabytes that a command line astap_index_solve on this machine is sharing,
rem and that cost minutes to rebuild. To take that as well:
rem
rem   uninstall.cmd -RemoveCache

setlocal
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0install.ps1" -Uninstall %*
set EXITCODE=%ERRORLEVEL%

rem Double-clicked means no arguments and a console of our own that is about to
rem close with the output still in it. Anything with arguments came from a shell
rem that will still be there afterwards, so never block that one.
if "%~1"=="" (
    echo %cmdcmdline% | find /i "%~nx0" >nul
    if not errorlevel 1 pause
)

exit /b %EXITCODE%
