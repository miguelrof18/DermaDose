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

## Cómo calcula

El lector obtiene del RDSR el kerma de referencia, la geometría, la colimación y la calidad del haz por evento. El motor sitúa los haces y la superficie en un sistema común, busca la primera intersección con la piel, aplica las correcciones dosimétricas configuradas y acumula los eventos. El visor muestra el mapa acumulado y la ubicación del máximo. La colocación estimada mediante TC y RDSR depende de las hipótesis de anclaje y mesa.

## Estructura

- `skindose_verifier/`: lector RDSR, reconstrucción, dosimetría, visor y datos requeridos.
- `slicer/DermaDoseSurface/`: módulo de extracción y exportación.
- `examples/`: instrucciones de demostración.
- `docs/`: configuración y colaboración.
- `tests/`: regresiones de lectura, geometría y cálculo.

Para contribuir, crea una rama y abre un PR. Ejecuta `python -m pip install -e ".[test]"` y `python -m pytest`. Incluye en el PR el caso reproducible y el efecto esperado. Usa casos sintéticos en GitHub y conserva los estudios de pacientes en el almacenamiento autorizado por el centro.
