# Actualización MRC 200 0407 (30 octubre 2026)

Se parte de la versión consolidada con pruebas físicas y caso reproducible. Los perfiles existentes y sus tablas se conservan. Esta entrega añade perfiles nominales explícitos, no reemplaza automáticamente configuraciones de ejemplos ni un ajuste basado en medidas.

## Uso

Desde la carpeta del programa:

```powershell
python -m pip install -e ".[test,spectral]"
python -m pytest
python reproducible/check.py
python -m skindose_verifier run RUTA_AL_RDSR --help
```

Los dos JSON nuevos están en `config/equipo_mrc200_conjunto_investigacion.json` (2.5 mm Al equivalente nominal) y `config/equipo_mrc200_conjunto_colimador_investigacion.json` (escenario aproximado 2.6 mm). Deben adaptarse la malla, el registro anatómico y geometría, así como la calibración y transmisión de mesa. No copiar estos números a ciegas al ejemplo07. El segundo perfil requiere verificar el limitador instalado. El primero no representa necesariamente la filtración total hasta el paciente.

`inherent_al_mm` se mantiene por compatibilidad del formato y representa aquí Al fijo efectivo del modelo. El nuevo bloque `tube_characterization` identifica claramente su significado. Los 1.7 mm Al del conjunto están incluidos en 2.5; no sumarlos otra vez. El código exige tablas de 11 grados y rechaza eventos >125 kV. La geometría ausente produce error en estos perfiles. Los filtros adicionales declarados en RDSR se conservan; el programa no inventa ni duplica 1 mm Al.

La distribución de los escenarios generados incluye Cu intermedios para interpolación; ello no implica que el equipo tenga esos filtros físicos. Los perfiles históricos 3.1 mm/12 grados y el ajuste condicionado previo permanecen disponibles para reproducir resultados antiguos. No se ha verificado la orientación anódica ni se ha implementado efecto talón o penumbra basada únicamente en el tamaño focal.

## Fuentes

- Philips, Service Information MRC 200 0407 ROT-GS1004, 4512-988-00132 REV AC (04/2011): sección4.4, Z-1, Z-2. Copia pública consultada: https://www.scribd.com/document/896599296/Tube . No se consiguió descargar el PDF original de servicio; no se distribuye una reconstrucción de él.
- Philips Azurion IFU Release2.2, 4523 001 01511: etiquetas p.416, filtros p.403. PDF oficial incluido en esta carpeta; URL https://www.documents.philips.com/assets/Instruction%20for%20Use/20210730/853843a4c24546639f5aad7500d34a2a.pdf?feed=ifu_docs_feed . Las etiquetas son ejemplos documentales y no fotografías de la sala.
- Los datos monoenergéticos Benmakhlouf y la procedencia de SpekPy se conservan en sus metadatos originales.

## Reproducción y limitaciones

Generación usada:

```powershell
python -m skindose_verifier.spectral build --out NUEVA_CARPETA --kvps 46,50,60,70,80,90,100,110,120,125 --al 2.5,2.6,3.5,3.6 --cu 0,0.1,0.2,0.3,0.4,0.6,0.9 --fields 5,10,20,25,35 --anode-angle 11
python scripts/compare_mrc200.py
```

El identificador y notas de procedencia del dataset distribuido se especializan como MRC200 nominal después de su generación; la tabla numérica no cambia. El CSV de sensibilidad compara a igual kVp/Cu/campo con 3.1 mm y 12 grados. En esa comparación se declara explícitamente 1 mm Al acompañante cuando Cu>0; no se deduce automáticamente en el lector RDSR. No son intervalos de incertidumbre ni validación clínica. Las HVL guardadas son calculadas, no medidas. Un mínimo de catálogo no es una medida espectral ni una garantía de exactitud.
