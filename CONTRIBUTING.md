# Trabajo con GitHub

La rama `main` contiene la versión integrada. Antes de un cambio, actualiza tu copia y crea una rama:

```powershell
git switch main
git pull --ff-only
git switch -c mejora/nombre-del-cambio
```

Después de editar, ejecuta las pruebas y revisa los archivos que se van a publicar:

```powershell
.\.venv\Scripts\python.exe -m pytest
git status
git diff
git add ruta/de/los/archivos
git commit -m "Descripción del cambio"
git push -u origin mejora/nombre-del-cambio
```

Abre un pull request hacia `main`. Revisa el cambio y las comprobaciones automáticas antes de integrarlo. Para sincronizar después, vuelve a `main` y ejecuta `git pull --ff-only`.

No se versionan entornos virtuales, resultados, copias de seguridad ni RDSR de estudios. Las muestras sintéticas se generan mediante el programa. Mantén el contenido binario y los hashes de los datasets científicos; no actualices un hash para ocultar un cambio de datos.

El historial anterior a esta incorporación se conserva en la etiqueta `respaldo-antes-importacion-20261012`. La memoria histórica continúa en `docs/`; la documentación de uso local actual está en `documentacion/`.
