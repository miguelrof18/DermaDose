# Guía sencilla para trabajar con DermaDose en GitHub

Esta guía describe el flujo recomendado sin necesidad de editar GitHub desde la línea de comandos. GitHub guarda el código y la documentación; cada ejecución se hace en tu ordenador.

## Primera instalación

1. Descarga el repositorio desde GitHub como ZIP o clónalo con GitHub Desktop.
2. Abre la carpeta `DermaDose` en VS Code.
3. Abre un terminal en esa carpeta y ejecuta:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
```

## Uso diario

1. Actualiza la copia local desde GitHub (`Pull` en GitHub Desktop o VS Code).
2. Trabaja en una rama con un nombre descriptivo, por ejemplo `mejora-lector-rdsr`.
3. Ejecuta `EJECUTAR_DEMO.bat` o `EJECUTAR_INVESTIGACION.bat`.
4. Comprueba los cambios con `python -m pytest`.
5. Guarda solo código, configuración, tests y documentación. No guardes `.venv`, resultados, RDSR de pacientes ni archivos generados.
6. Haz commit con un mensaje claro y publica la rama.
7. Abre un pull request hacia `main`. Explica qué cambias, cómo lo probaste y cualquier limitación.

## Qué debe ir y qué no debe ir al repositorio

Sí debe ir: `skindose_verifier/`, `config/`, `examples/`, `scripts/`, `tests/`, `documentacion/`, `docs/`, los lanzadores y los ficheros de configuración del repositorio.

No debe ir: `.venv/`, `resultados/`, `datos/` local, `archivo_legacy/`, DICOM/RDSR reales, mapas HTML generados, PNG, logs, claves ni secretos. El `.gitignore` ya los excluye.

## Para aceptar un cambio

Revisa el diff y la pestaña `Checks` del pull request. La integración continua ejecuta las pruebas en Windows y Ubuntu. Solo después de que las comprobaciones estén en verde conviene pulsar `Merge pull request`. La etiqueta `respaldo-antes-importacion-20261012` permite recuperar el punto anterior a la importación.

## Si algo sale mal

No borres la carpeta de trabajo. Vuelve a la rama `main`, crea otra rama desde el último estado conocido y conserva el pull request anterior para poder comparar. Para un error reproducible, usa la plantilla de incidencia y no adjuntes RDSR con datos identificativos.
