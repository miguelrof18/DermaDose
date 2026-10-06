# DermaDose

DermaDose es un prototipo de investigacion en Python para estimar y verificar mapas de dosis de entrada en piel en procedimientos intervencionistas a partir de eventos DICOM RDSR. El programa reconstruye la geometria del arco, proyecta el campo de radiacion sobre una superficie de paciente o fantoma, acumula contribuciones por evento y genera resultados trazables con mapas 2D, visualizacion 3D y archivos JSON/NPZ.

El codigo procede de `skindose-verifier` v0.5.3 y se ha reorganizado para trabajar como repositorio GitHub mantenible.

## Estado del proyecto

Este software esta en fase de investigacion. Las configuraciones de demostración son sintéticas y no representan una calibración clínica; las tablas publicadas se identifican por separado con su procedencia y dominio. No debe utilizarse para tomar decisiones clinicas sin validacion fisica, caracterizacion local del equipo y revision por radiofisica.

## Estructura

```text
DermaDose/
  src/skindose_verifier/      Codigo fuente del paquete Python
  examples/                   Plantillas y entradas sinteticas
  tests/                      Pruebas automaticas
  scripts/                    Utilidades de reproduccion y medicion
  scripts/windows/            Lanzadores .bat heredados
  docs/                       Memoria, auditoria tecnica, validacion y roadmap
  pyproject.toml              Configuracion del paquete instalable
  requirements.txt            Dependencias principales
```

Los resultados pesados generados en versiones previas no se han incluido en el commit inicial para mantener el repositorio limpio. Pueden regenerarse con los comandos del paquete.

## Instalacion rapida

En Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
.\.venv\Scripts\python.exe -m pytest
.\.venv\Scripts\python.exe -m skindose_verifier demo --out demo_output --stage skin_dose
```

En Linux o macOS:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e ".[test]"
.venv/bin/python -m pytest
.venv/bin/python -m skindose_verifier demo --out demo_output --stage skin_dose
```

Tambien queda instalado el comando `dermadose`, equivalente al comando historico `skindose`.

## Comandos utiles

```bash
python -m skindose_verifier inspect archivo_anonimizado.dcm
python -m skindose_verifier demo --out demo_elipse --model ellipse --stage skin_dose
python -m skindose_verifier demo --out demo_mesh --model mesh --stage skin_dose
python -m skindose_verifier run archivo_anonimizado.dcm --config equipo_validado.json --out salida --stage skin_dose
python -m skindose_verifier convergence archivo_anonimizado.dcm --config equipo_validado.json --stage skin_dose --resolutions 10,5,2.5
```

La documentacion tecnica completa esta en `docs/README_original_v0.5.3.md`, junto con la memoria cientifica y los documentos de validacion y roadmap.

## Licencias y procedencia de assets

Las mallas de demostracion incluidas en `src/skindose_verifier/data/` conservan sus notas de licencia en:

- `src/skindose_verifier/data/XCAT_LICENSE.txt`
- `src/skindose_verifier/data/PYSKINDOSE_LICENSE.txt`
- `docs/ASSET_ORIGEN.md`

## Proximos pasos tecnicos

Ya se incluyen las tablas del informe KS-ASF-201101-IR (Benmakhlouf et al., 2011), con interpolación acotada y procedencia por valor. El proveedor es opcional y produce dosis superficial equivalente en **agua**; la demostración histórica mantiene sus factores sintéticos.

- [Informe técnico, límites físicos y resultados de verificación](docs/BENMAKHLOUF_2011.md).
- [Configuración de investigación con factores publicados](examples/equipo_benmakhlouf_investigacion.json).
- [Datos originales y metadatos](src/skindose_verifier/data/benmakhlouf_2011/).

El proveedor rechaza extrapolaciones y combinaciones con `bsf`/`f_skin` heredados. Los factores de mesa, calibración y geometría del ejemplo siguen siendo sintéticos. La incorporación de tablas publicadas no constituye validación clínica.
