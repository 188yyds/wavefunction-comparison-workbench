@echo off
setlocal
cd /d "%~dp0"
if exist "%~dp0.venv\Scripts\pythonw.exe" (
  start "" "%~dp0.venv\Scripts\pythonw.exe" "%~dp0launcher.pyw"
  exit /b 0
)
where pythonw >nul 2>nul
if not errorlevel 1 (
  start "" pythonw "%~dp0launcher.pyw"
  exit /b 0
)
where pyw >nul 2>nul
if not errorlevel 1 (
  start "" pyw -3 "%~dp0launcher.pyw"
  exit /b 0
)
where python >nul 2>nul
if not errorlevel 1 (
  start "" python "%~dp0launcher.pyw"
  exit /b 0
)
echo Python with PyQt5, Matplotlib and h5py is required. See README.md.
pause
