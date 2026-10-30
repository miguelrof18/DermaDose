@echo off
setlocal
cd /d "%~dp0"
if exist "resultados\research_xcat\mapa_3d_interactivo.html" (
  start "" "resultados\research_xcat\mapa_3d_interactivo.html"
  exit /b 0
)
if exist "resultados\demo_mesh\mapa_3d_interactivo.html" (
  start "" "resultados\demo_mesh\mapa_3d_interactivo.html"
  exit /b 0
)
echo Ejecuta primero EJECUTAR_DEMO.bat o EJECUTAR_INVESTIGACION.bat.
pause
exit /b 1
