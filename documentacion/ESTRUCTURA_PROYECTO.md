# Estructura del proyecto

## Zonas que se deben editar

- `skindose_verifier/`: código Python del cálculo, lector RDSR, geometría, ray tracing, factores, métricas e informes.
- `config/`: configuraciones reutilizables del equipo y del fantoma.
- `tests/`: pruebas automáticas. Deben ejecutarse después de cualquier cambio del núcleo.
- `scripts/`: comprobaciones, validación y reproducción de estudios.
- `documentacion/`: guías de uso y documentación técnica.

## Zonas de datos y resultados

- `datos/`: datasets versionados, mallas, tablas y metadatos SHA-256. No editar manualmente las tablas protegidas.
- `resultados/`: salidas generadas por demos, investigaciones y pruebas. No es código fuente.
- `examples/`: configuraciones y entradas sintéticas de referencia.
- `archivo_legacy/`: respaldos y material histórico; no se utiliza en la ejecución normal.

## Punto de entrada para cada tarea

- `EJECUTAR_DEMO.bat`: demostración sintética.
- `EJECUTAR_INVESTIGACION.bat`: ejecución de investigación.
- `EJECUTAR_EJEMPLO.bat`: lee las carpetas externas `..\ejemplos rdsr\ejemploXX` y guarda allí los resultados.

Los lanzadores de la raíz se conservan por compatibilidad con el Explorador de Windows. Las copias de `scripts/windows/` son auxiliares para usuarios técnicos.

## Regla de mantenimiento

Una configuración nueva debe colocarse en `config/`, un dataset en `datos/` con su metadato, y una salida en una subcarpeta de `resultados/`. No se deben guardar DICOM de prueba ni resultados dentro de `skindose_verifier/`.

## Comprobación antes de entregar cambios

```powershell
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m skindose_verifier demo --model mesh --out resultados\demo_mesh --stage skin_dose
```

