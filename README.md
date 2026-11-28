# DermaDose

DermaDose estima la distribución de dosis de entrada sobre una superficie cutánea y su máximo (PSD) en procedimientos intervencionistas. Es un programa de investigación sin validación clínica.

## Empezar

Requiere Python 3.10 o posterior. En PowerShell, desde la raíz del repositorio:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\python.exe -m skindose_verifier demo --model mesh --out resultados/demo --stage skin_dose
```

Abre `resultados/demo/mapa_3d_interactivo.html` en el navegador y consulta `resultados/demo/resultados.json` para ver la PSD y la trazabilidad.

## Ejecutar un caso

Prepara una carpeta con un RDSR DICOM y un solo JSON de equipo, copiado de [los perfiles disponibles](config/). Para utilizar la superficie del paciente, añade un STL y su JSON homónimo exportados por DermaDose Surface.

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier.example_runner "C:\casos\ejemplo" --mesh patient
```

Los resultados se guardan en la subcarpeta `resultados` del caso. Consulta [el ejemplo](examples/ejemplo/README.md) y la [guía de equipo](docs/configuracion.md).

## DermaDose Surface para 3D Slicer

En Slicer, añade `slicer/DermaDoseSurface` a las rutas de módulos adicionales en **Edit → Application Settings → Modules** y reinicia. Abre **DermaDose → DermaDose Surface**, selecciona el volumen de TC, extrae y revisa el contorno y exporta el STL con su JSON. Conserva ambos juntos. La superficie exportada describe el contorno externo del TC.

## Colocación del paciente

La malla se sitúa con una sola cadena: TC (LPS) → tablero (origen en el borde de cabecera) → isocentro, con las lecturas de mesa de cada evento. Los eslabones que no se han medido (registro del TC en el tablero, referencia tablero–isocentro) los estima el lanzador con hipótesis explícitas que quedan en los resultados. Consulta [Sistemas de referencia](docs/sistemas_de_referencia.md) y la herramienta `python -m skindose_verifier.table_calibration` para caracterizar los movimientos de la mesa.

## Verificación

La [verificación RFM 2024](docs/verificacion_rfm2024.md) compara DermaDose con las 27 medidas del conjunto de Gañán Mora et al. (Rev Fis Med 2024). Todos los casos quedan dentro de la tolerancia del 13 %.

## Cómo calcula

El lector obtiene del RDSR el kerma de referencia, la geometría, la colimación y la calidad del haz por evento. El motor sitúa los haces y la superficie en un sistema común, busca la primera intersección con la piel, aplica las correcciones dosimétricas configuradas y acumula los eventos. El visor muestra el mapa acumulado y la ubicación del máximo. La colocación estimada mediante TC y RDSR depende de las hipótesis de anclaje y mesa.

## Qué contiene el resultado

`resultados.json` identifica la magnitud por su etapa: `max_air_kerma_gy` en `free_air` y `calibrated_air`, y `max_modeled_surface_dose_gy` en `skin_dose`, junto con `psd_point_definition` y `dose_medium`. El campo `psd_point_gy` se mantiene por compatibilidad con valor idéntico. La salida en agua no es dosis a una profundidad cutánea concreta. El bloque `performance` desglosa el tiempo por fases. `configuration_warnings` avisa de claves del JSON que el programa no utiliza.

## Estructura

- `skindose_verifier/`: lector RDSR, reconstrucción, dosimetría, visor y datos requeridos.
- `slicer/DermaDoseSurface/`: módulo de extracción y exportación.
- `examples/`: instrucciones de demostración.
- `docs/`: configuración y colaboración.
- `tests/`: regresiones de lectura, geometría y cálculo.

Para contribuir, crea una rama y abre un PR. Ejecuta `python -m pip install -e ".[test]"` y `python -m pytest`. Los cambios de cada versión están en [CHANGELOG.md](CHANGELOG.md). Incluye en el PR el caso reproducible y el efecto esperado. Usa casos sintéticos en GitHub y conserva los estudios de pacientes en el almacenamiento autorizado por el centro.
