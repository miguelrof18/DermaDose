@echo off
setlocal
cd /d "%~dp0" || goto error
set "PY=.venv\Scripts\python.exe"
if not exist "%PY%" (
    python -m venv .venv || goto error
)
"%PY%" -m pip show dermadose >nul 2>&1
if errorlevel 1 (
    "%PY%" -m pip install -e . || goto error
)
set "CASE=%~dp0examples\ejemplo"
echo Introduce otra ruta o pulsa Intro para usar el ejemplo incluido.
echo Pulsa Intro para usar: "%CASE%"
set /p "CASE=Carpeta: "
set "CASE=%CASE:"=%"
if not exist "%CASE%\." goto error
echo.
"%PY%" -m skindose_verifier.example_runner "%CASE%" --interactive --open-viewer || goto error
echo.
pause
exit /b 0

:error
echo.
echo No se pudo completar la ejecucion. Revisa el mensaje anterior.
pause
exit /b 1
