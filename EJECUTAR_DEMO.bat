@echo off
setlocal
cd /d "%~dp0"
if errorlevel 1 goto error

rem Reutiliza un entorno ya creado, o encuentra un Python sin depender de py.
if exist ".venv\Scripts\python.exe" goto run
set "PY_CMD="
python -c "import sys; sys.exit(sys.version_info < (3, 10))" >nul 2>&1
if not errorlevel 1 set "PY_CMD=python"
if not defined PY_CMD (
    py -3 -c "import sys; sys.exit(sys.version_info < (3, 10))" >nul 2>&1
    if not errorlevel 1 set "PY_CMD=py -3"
)
if not defined PY_CMD (
    python3 -c "import sys; sys.exit(sys.version_info < (3, 10))" >nul 2>&1
    if not errorlevel 1 set "PY_CMD=python3"
)
if not defined PY_CMD goto no_python

echo Creando el entorno virtual con %PY_CMD%...
%PY_CMD% -m venv .venv
if errorlevel 1 goto error

:run
".venv\Scripts\python.exe" -c "import sys; sys.exit(sys.version_info < (3, 10))" >nul 2>&1
if errorlevel 1 goto bad_venv
".venv\Scripts\python.exe" -m pip install -e ".[test]"
if errorlevel 1 goto error
".venv\Scripts\python.exe" -m pytest
if errorlevel 1 goto error
".venv\Scripts\python.exe" -m skindose_verifier demo --out resultados\demo_output --stage skin_dose
if errorlevel 1 goto error
".venv\Scripts\python.exe" -m skindose_verifier demo --model mesh --mesh-edge-mm 10 --out resultados\demo_mesh --stage skin_dose
if errorlevel 1 goto error
echo.
echo LISTO: abre resultados\demo_mesh\mapa_3d_interactivo.html en el navegador
pause
exit /b 0
:no_python
echo.
echo No se ha encontrado Python 3.10 o posterior.
echo Instala Python desde https://www.python.org/downloads/
echo Cierra VS Code por completo y vuelve a abrirlo despues de instalarlo.
echo Luego haz doble clic de nuevo en EJECUTAR_DEMO.bat.
pause
exit /b 1
:bad_venv
echo.
echo El Python de .venv no funciona o es anterior a 3.10.
echo Borra la carpeta .venv y ejecuta otra vez EJECUTAR_DEMO.bat.
pause
exit /b 1
:error
echo.
echo Hubo un error. Copia el texto de esta ventana para diagnosticarlo.
pause
exit /b 1
