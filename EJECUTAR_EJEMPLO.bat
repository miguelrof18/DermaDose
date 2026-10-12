@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo No se encuentra el entorno virtual. Ejecuta primero EJECUTAR_DEMO.bat.
  pause
  exit /b 1
)
set /p NUM=Escribe el numero del ejemplo (por ejemplo 01):
rem Los ejemplos de entrada se mantienen fuera del programa, en la carpeta hermana.
set "INPUT_DIR=..\ejemplos rdsr\ejemplo%NUM%"
set "OUTPUT_DIR=%INPUT_DIR%\resultados"
if not exist "%INPUT_DIR%" (
  echo No existe la carpeta externa de ejemplo: %INPUT_DIR%
  pause
  exit /b 1
)
set "RDSR_FILE="
for /f "delims=" %%F in ('dir /b /a-d "%INPUT_DIR%\*.dcm" 2^>nul') do if not defined RDSR_FILE set "RDSR_FILE=%INPUT_DIR%\%%F"
if not defined RDSR_FILE (
  echo No se encontro ningun archivo DICOM .dcm en %INPUT_DIR%
  pause
  exit /b 1
)
set "CONFIG_FILE="
for /f "delims=" %%F in ('dir /b /a-d "%INPUT_DIR%\*.json" 2^>nul') do if not defined CONFIG_FILE set "CONFIG_FILE=%INPUT_DIR%\%%F"
if not defined CONFIG_FILE (
  echo No se encontro ningun archivo de configuracion .json en %INPUT_DIR%
  pause
  exit /b 1
)
echo Ejecutando ejemplo%NUM%...
echo RDSR: %RDSR_FILE%
echo Configuracion: %CONFIG_FILE%
".venv\Scripts\python.exe" -m skindose_verifier run "%RDSR_FILE%" --config "%CONFIG_FILE%" --out "%OUTPUT_DIR%" --stage skin_dose
if errorlevel 1 (
  echo El ejemplo no se pudo ejecutar.
  pause
  exit /b 1
)
echo Resultado guardado en %OUTPUT_DIR%
start "" "%OUTPUT_DIR%\mapa_3d_interactivo.html"
pause
