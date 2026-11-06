@echo off
cd /d "%~dp0"
title Fetch reference data 2
set "PY=%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
if not exist "%PY%" set "PY=python"
"%PY%" -I "%~dp0fetch_refs2.py"
echo.
pause
