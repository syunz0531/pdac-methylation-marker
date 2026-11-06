@echo off
cd /d "%~dp0"
title PDAC extra extraction
set "PY=%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
if not exist "%PY%" set "PY=python"
"%PY%" "%~dp0extract_extra.py"
echo.
pause
