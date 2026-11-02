@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo No se encuentra el entorno virtual. Ejecuta primero EJECUTAR_DEMO.bat.
  pause
  exit /b 1
)
".venv\Scripts\python.exe" -m skindose_verifier.example_runner --examples-root "..\ejemplos rdsr" --interactive --open-viewer
if errorlevel 1 (
  echo El ejemplo no se pudo ejecutar. Revisa el mensaje anterior.
  pause
  exit /b 1
)
echo Resultado guardado en la carpeta resultados del ejemplo elegido.
pause
