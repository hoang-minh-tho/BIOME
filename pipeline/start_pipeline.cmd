@echo off
cd /d "%~dp0"
where py >nul 2>nul
if errorlevel 1 (
    python "%~dp0run_pipeline.py"
) else (
    py -3 "%~dp0run_pipeline.py"
)
if errorlevel 1 (
    echo.
    echo Install Python 3.9 or newer with Tcl/Tk from https://www.python.org/downloads/windows/
    echo Then run this file again. See HUONG_DAN.md for instructions.
    pause
)
