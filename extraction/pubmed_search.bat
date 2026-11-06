@echo off
cd /d "%~dp0"
title PubMed search
set "PY=%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
if not exist "%PY%" set "PY=python"
"%PY%" -I "%~dp0pubmed_search.py"
echo.
pause
