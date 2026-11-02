# DermaDose

Calculador de mapas de dosis superficial a partir de informes DICOM RDSR. Versión local incorporada en octubre de 2026: XCAT50, factores históricos o Benmakhlouf, lector ampliado, visor 3D con isodosis y ray tracing optimizado.

## Instalar y ejecutar

Necesitas Python 3.10 o posterior. Las comprobaciones automáticas usan Python 3.11.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
.\.venv\Scripts\python.exe -m skindose_verifier demo --model mesh --out resultados\demo_mesh --stage skin_dose
```

En Windows puedes usar los lanzadores con doble clic:
- `EJECUTAR_DEMO.bat`: prepara el entorno, ejecuta las pruebas y genera demos.
- `EJECUTAR_INVESTIGACION.bat`: reproduce el estudio completo (tarda más).
- `VER_MAPA_XCAT.bat`: abre el mapa ya generado de investigación o demo.
- `EJECUTAR_EJEMPLO.bat`: ejecuta un ejemplo externo.

En Linux/macOS usa `python3 -m venv .venv` y `.venv/bin/python` en los comandos equivalentes.

## Ejemplos externos

El lanzador busca una carpeta hermana del repositorio:

```text
Calculador_Dosis_Piel_Proyecto_v0.5.3/
  DermaDose/
  ejemplos rdsr/
    ejemplo01/
      nombre_del_equipo.dcm
      config.json
      resultados/
```

El lanzador pregunta primero el ejemplo y después si usar XCAT o el STL del paciente guardado en esa carpeta, con cualquier nombre. Si hay varios candidatos, permite elegir. Solo se necesitan RDSR, JSON de equipo y, para paciente, STL. La colocación provisional del equipo se muestra en texto y en una ficha visual sin pedir confirmación ni `paciente.json`. Véase [ejemplos con malla local](documentacion/EJEMPLOS_CON_MALLA_LOCAL.md). Las salidas quedan en `ejemploXX/resultados/` y se actualizan al repetir.

Para seleccionar archivos explícitamente:

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier run "ruta\estudio.dcm" --config "ruta\equipo.json" --out "ruta\resultados" --stage skin_dose
```

Las rutas del JSON se interpretan desde su carpeta. Al copiarlo junto al RDSR, adapta `patient.mesh_path`. Desde `ejemplos rdsr/ejemplo01/` sería `../../DermaDose/skindose_verifier/data/xcat_wholebody_body.stl` para esta copia del programa.

## Configuraciones

- `config/equipo_xcat_sintetico_historico.json`: XCAT50 y factores constantes sintéticos.
- `config/equipo_xcat_benmakhlouf_investigacion.json`: XCAT50 y proveedor Benmakhlouf.
- [Guía de campos JSON](documentacion/GUIA_CONFIGURACION_JSON.md). Para la política concreta de compatibilidad/fallback prevalecen el código y la auditoría de la ejecución.

Las configuraciones incluyen geometría y calibraciones sintéticas. El software es un prototipo de investigación sin validación clínica.

## Estructura

```text
skindose_verifier/  Código principal
  data/            Mallas, licencias, plantilla HTML y tablas SHA-256
config/            Configuraciones reutilizables
examples/          Plantillas de referencia
scripts/           Validación y estudios
tests/             Pruebas automáticas
documentacion/     Documentación de la versión local
docs/              Memoria y documentación anterior conservadas
.github/workflows/ Comprobaciones Windows y Linux
```

Los entornos, DICOM de estudios, resultados, archivos ICRP, respaldos y datos voluminosos permanecen locales. Los recursos necesarios y redistribuibles de la antigua carpeta `datos/` están en `skindose_verifier/data/`.

## Desarrollo

```powershell
.\.venv\Scripts\python.exe -m pytest
.\.venv\Scripts\python.exe scripts\validate_benmakhlouf.py
```

Consulta [CONTRIBUTING.md](CONTRIBUTING.md) para trabajar con ramas y pull requests. La etiqueta `respaldo-antes-importacion-20261012` conserva la versión anterior de GitHub. El paquete mantiene la versión 0.5.3 y cada commit identifica los cambios.

## Procedencia

Se conservan [licencia XCAT](skindose_verifier/data/XCAT_LICENSE.txt), [licencia PySkinDose](skindose_verifier/data/PYSKINDOSE_LICENSE.txt) y [procedencia de mallas](documentacion/ASSET_ORIGEN.md). Las tablas conservan los bytes y hashes de la copia local verificada. No se añade una licencia general nueva al código.
