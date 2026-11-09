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
"%PY%" -m skindose_verifier demo --model mesh --out resultados/demo --stage skin_dose || goto error
start "" "resultados\demo\mapa_3d_interactivo.html"
echo.
pause
exit /b 0

:error
echo.
echo No se pudo completar la ejecucion. Revisa el mensaje anterior.
pause
exit /b 1
