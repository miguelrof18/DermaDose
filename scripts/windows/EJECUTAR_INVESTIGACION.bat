@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
    echo Primero ejecuta EJECUTAR_DEMO.bat para instalar el entorno.
    pause
    exit /b 1
)
".venv\Scripts\python.exe" -m pip install -e ".[test]"
if errorlevel 1 goto error
".venv\Scripts\python.exe" scripts\reproduce_research.py
if errorlevel 1 goto error
echo Estudio terminado: abre research_results y results\research_xcat.
start "" "results\research_xcat\mapa_3d_interactivo.html"
pause
exit /b 0
:error
echo El estudio se ha detenido. Copia el error que aparece encima.
pause
exit /b 1
