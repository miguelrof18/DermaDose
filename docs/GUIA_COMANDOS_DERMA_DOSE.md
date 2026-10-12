# Guía práctica de DermaDose y GitHub
Guía diaria para Windows y PowerShell.
## Preparar el entorno
```powershell
cd "C:\Users\Miguel\OneDrive\Escritorio\Calculador_Dosis_Piel_Proyecto_v0.5.3\DermaDose"
git status
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
```
En VS Code selecciona `.venv\Scripts\python.exe` como intérprete.
## Ejecutar DermaDose
```powershell
.\EJECUTAR_DEMO.bat
.\EJECUTAR_INVESTIGACION.bat
.\EJECUTAR_EJEMPLO.bat
```
El último lanzador pide `01`, `02`, etc. y busca el primer `.dcm` y `.json` en `..\ejemplos rdsr\ejemploXX\`. Los resultados quedan en la carpeta `resultados` del ejemplo.
## Actualizar desde GitHub
```powershell
git switch main
git pull --ff-only
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
```
## Trabajar en una mejora
```powershell
git switch main
git pull --ff-only
git switch -c mejora-nombre-descriptivo
```
Después de editar y probar:
```powershell
git status
git diff
git add .
git commit -m "Describe brevemente el cambio"
git push -u origin mejora-nombre-descriptivo
```
En GitHub pulsa **Compare & pull request** y crea el pull request hacia `main`. Revisa los **Checks** antes de pulsar **Merge pull request**.
## Pruebas
```powershell
 .\.venv\Scripts\python.exe -m pytest
 .\.venv\Scripts\python.exe scripts\validate_benmakhlouf.py
```
## Ramas y recuperación
```powershell
git branch
git branch -a
git log --oneline --decorate -10
git switch main
```
`git restore ruta\archivo` descarta cambios de un archivo concreto. No uses `git reset --hard` sin confirmar antes qué se perderá.
## Qué se sube
Se suben código, configuraciones, tests, scripts y documentación. No se suben `.venv`, `resultados`, mapas, logs, RDSR privados, claves ni respaldos; `.gitignore` ya los excluye.
