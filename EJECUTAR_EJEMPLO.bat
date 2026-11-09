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
echo 1. Factores sinteticos historicos
echo 2. Tablas publicadas de Benmakhlouf
echo 3. Philips Azurion: Benmakhlouf y SpekPy
choice /c 123 /n /m "Elige el perfil de calculo: "
if errorlevel 4 goto error
if not errorlevel 1 goto error
if errorlevel 1 set "CONFIG=%~dp0skindose_verifier\data\equipo_sintetico.json"
if errorlevel 2 set "CONFIG=%~dp0examples\equipo_benmakhlouf_investigacion.json"
if errorlevel 3 set "CONFIG=%~dp0examples\ejemplo\equipo_azurion.json"
"%PY%" -m skindose_verifier.example_runner "%CASE%" --config "%CONFIG%" --interactive --open-viewer || goto error
echo.
pause
exit /b 0

:error
echo.
echo No se pudo completar la ejecucion. Revisa el mensaje anterior.
pause
exit /b 1
