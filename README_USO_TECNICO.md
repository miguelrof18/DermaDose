# Uso técnico

El paquete ejecutable es `skindose_verifier`. El entorno recomendado es `.venv`.

## Comandos principales

```powershell
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m skindose_verifier demo --model mesh --mesh-edge-mm 10 --out resultados\demo_manual --stage skin_dose
.\.venv\Scripts\python.exe -m skindose_verifier run resultados\research_xcat\rdsr_sintetico.dcm --config config\equipo_xcat_benmakhlouf_investigacion.json --out resultados\xcat_benmakhlouf --stage skin_dose
```

`config` contiene configuraciones editables; `datos` es una copia navegable de los
datasets; la copia canónica usada por el paquete permanece en
`skindose_verifier/data`. Los resultados nuevos deben escribirse bajo `resultados`.

El dataset Benmakhlouf verifica SHA-256 antes de calcular y rechaza extrapolaciones.
## Fantomas comparables

XCAT50 sigue siendo el fantoma predeterminado de las configuraciones existentes. Una configuración puede declarar explícitamente un fantoma alternativo:

```json
"phantom": {
  "name": "ICRP_AF",
  "mesh_path": "../datos/fantomas/icrp_adult_female.stl",
  "mesh_unit": "mm",
  "mesh_axes_matrix": [[1,0,0],[0,1,0],[0,0,1]],
  "mesh_origin_mm": [0,0,0]
}
```

Los nombres admitidos son `XCAT50`, `ICRP_AM` e `ICRP_AF`. El programa exige una malla STL/OBJ/PLY para los dos perfiles ICRP; no incluye mallas ICRP con licencia en esta distribución. No se debe renombrar otra malla como ICRP.

Para ejecutar un mismo RDSR en XCAT50 y un perfil alternativo:

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier compare-phantoms `
  ruta\procedimiento.dcm `
  --xcat-config examples\equipo_xcat_benmakhlouf_investigacion.json `
  --alternative-config config\equipo_icrp_af.json `
  --stage skin_dose `
  --out resultados\comparacion_fantomas
```

El informe `comparacion_fantomas/comparacion_fantomas.json` incluye PSD, dosis media, D2/D5/D10, área irradiada, integral dosis-área, diferencia relativa de PSD, distancia entre superficies y diferencia espacial de dosis. La comparación usa el mismo RDSR y la misma cadena dosimétrica; las diferencias se atribuyen al muestreo y a la anatomía de las mallas, pero no constituyen validación clínica.
