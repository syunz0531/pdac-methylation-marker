@echo off
cd /d "%~dp0"
title PDAC methylation data extraction

set "PY="
set "PYUSER=%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
if exist "%PYUSER%" set "PY=%PYUSER%"
if not defined PY (
  python -c "import sys; sys.exit(0 if sys.version_info>=(3,8) else 1)" >nul 2>nul && set "PY=python"
)

if not defined PY (
  echo [1/2] Python not found. Installing Python 3.12 with winget...
  echo       If asked to agree to terms, type Y and press Enter.
  winget install -e --id Python.Python.3.12 --scope user
  if exist "%PYUSER%" set "PY=%PYUSER%"
)

if not defined PY (
  echo.
  echo Could not find Python after install. Please show this window to Claude.
  pause
  exit /b 1
)

echo [2/2] Downloading and extracting data. This can take hours. Keep this window open.
"%PY%" -I "%~dp0extract.py"
echo.
pause
