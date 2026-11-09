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
set "CASE=%USERPROFILE%\OneDrive\Escritorio\ejemplos rdsr\ejemplo08"
echo Introduce la ruta de la carpeta con el RDSR.
echo Pulsa Intro para usar: "%CASE%"
set /p "CASE=Carpeta: "
set "CASE=%CASE:"=%"
if not exist "%CASE%\." goto error
choice /c 12 /n /m "Superficie: 1 = XCAT, 2 = TC del paciente: "
if errorlevel 3 goto error
set "MESH=xcat"
if errorlevel 2 set "MESH=patient"
"%PY%" -m skindose_verifier.example_runner "%CASE%" --config "%CD%\examples\ejemplo\equipo_azurion.json" --mesh %MESH% --open-viewer || goto error
echo.
pause
exit /b 0

:error
echo.
echo No se pudo completar la ejecucion. Revisa el mensaje anterior.
pause
exit /b 1
