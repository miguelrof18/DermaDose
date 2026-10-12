# Guía del archivo JSON de configuración

## 1. Qué archivos intervienen

Una ejecución utiliza dos entradas principales:

1. **RDSR DICOM** (`.dcm`): contiene los eventos registrados por el equipo: kerma o KAP, kVp, filtración, ángulos, SID/SOD/SRP, tamaño de campo, mesa y fabricante/modelo.
2. **Configuración JSON** (`.json`): describe cómo interpretar esos eventos y cómo reconstruir la geometría y la dosis sobre el paciente o fantoma.

El JSON no sustituye al RDSR ni debe editarse para corregir datos medidos. Si un campo no está en el RDSR, solo puede completarse mediante una hipótesis documentada en la configuración o mediante una caracterización del equipo.

Ejemplo de ejecución:

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier run `
  mi_estudio.dcm `
  --config config\mi_equipo.json `
  --out resultados\mi_estudio `
  --stage skin_dose
```

Las rutas relativas del JSON se resuelven respecto a la carpeta donde está el propio JSON.

## 2. Campos generales

```json
{
  "version": "MI-EQUIPO-v1",
  "manufacturer": "FABRICANTE",
  "model": "MODELO",
  "validation_status": "estado de validación",
  "dose_medium": "skin"
}
```

- `version`: identificador de la configuración. Debe cambiarse cuando cambie una calibración, geometría o dataset.
- `manufacturer` y `model`: deben coincidir con los valores del RDSR. El programa rechaza una combinación diferente para evitar aplicar la geometría de un equipo a otro.
- `validation_status`: texto de trazabilidad; por ejemplo, `synthetic_self_consistency_only` o `research_not_clinically_validated`.
- `dose_medium`: `skin` para factores históricos directos o `water` cuando se utiliza el proveedor tabulado de Benmakhlouf.

## 3. Bloque `geometry`

Este bloque define la geometría fija del equipo y de la mesa:

```json
"geometry": {
  "source_to_isocenter_mm": 700,
  "source_to_detector_mm": 1100,
  "source_to_irp_mm": 550,
  "field_reference": "detector",
  "primary_sign": 1,
  "secondary_sign": 1,
  "rotation_order": "primary_then_secondary",
  "isocenter_mm": [0, 0, 0],
  "table_coordinates": "isocenter",
  "table_origin_mm": [0, 0, 0],
  "table_matrix": [[1,0,0],[0,1,0],[0,0,1]],
  "table_plane_z_mm": -125,
  "table_half_width_mm": 300,
  "table_half_length_mm": 600
}
```

- `source_to_isocenter_mm`: distancia fuente-isocentro.
- `source_to_detector_mm`: SID. Debe ser mayor que la distancia fuente-isocentro.
- `source_to_irp_mm`: distancia de referencia usada por el equipo o por la calibración.
- `field_reference`: referencia del tamaño de campo: `detector` o `isocenter`.
- `primary_sign`, `secondary_sign`: signos de los ángulos del RDSR al transformarlos al sistema del programa.
- `rotation_order`: orden de las rotaciones, `primary_then_secondary` o `secondary_then_primary`.
- `isocenter_mm`: origen del sistema clínico.
- `table_coordinates`: sistema de coordenadas de la mesa, `isocenter` o `legacy`.
- `table_origin_mm`: origen de la mesa en el sistema elegido.
- `table_matrix`: matriz ortonormal que transforma los ejes de la mesa.
- `table_plane_z_mm`, `table_half_width_mm`, `table_half_length_mm`: plano y dimensiones de la mesa en milímetros.

Si hay geometrías diferentes por plano (por ejemplo, frontal y lateral), pueden declararse en `geometry.planes`. Cada plano debe estar caracterizado; no se recomienda inventar una geometría para ocultar datos faltantes.

## 4. Bloque `patient` y fantomas

Para una superficie plana de prueba:

```json
"patient": {
  "type": "plane",
  "width_mm": 360,
  "length_mm": 360,
  "z_mm": -100,
  "normal": [0, 0, -1],
  "thickness_mm": 220,
  "origin_on_table_mm": [0, 0, 0],
  "resolution_mm": 5
}
```

Para una elipse se usa `type: "ellipse"` con sus dimensiones geométricas. Para XCAT o una malla de superficie se usa:

```json
"patient": {
  "type": "mesh",
  "mesh_path": "../datos/fantomas/xcat_surface.stl",
  "mesh_unit": "mm",
  "mesh_axes_matrix": [[1,0,0],[0,0,1],[0,-1,0]],
  "mesh_origin_mm": [0, -500, 40],
  "origin_on_table_mm": [0, 0, 0],
  "resolution_mm": 5,
  "mesh_max_edge_mm": 10
}
```

También puede seleccionarse un fantoma con el bloque superior:

```json
"phantom": {
  "name": "XCAT50"
}
```

Los nombres admitidos por el programa son `XCAT50`, `ICRP_AM` e `ICRP_AF`; los dos últimos requieren `mesh_path`. La malla y su SHA-256 deben conservarse sin modificaciones cuando se pretenda reproducir un resultado.

## 5. Bloque `factors`

Cada factor debe documentar su procedencia e identificador. Puede ser constante:

```json
"cf": {
  "value": 1.0,
  "provenance": "calibración del equipo, fecha...",
  "identifier": "CF-v1",
  "relative_standard_uncertainty": 0.05
}
```

O puede ser una tabla multidimensional interpolada:

```json
"bsf": {
  "axes": {
    "kvp": [70, 80, 90],
    "field_area_cm2": [100, 400]
  },
  "grid": [[1.20, 1.25], [1.24, 1.30], [1.28, 1.35]],
  "provenance": "medición o publicación",
  "identifier": "BSF-v1"
}
```

El orden de `grid` debe corresponder al orden de los ejes. El programa interpola únicamente dentro del dominio tabulado y rechaza la extrapolación. Los factores positivos y los ejes estrictamente crecientes son obligatorios.

Los factores históricos habituales son:

- `cf`: conversión de KAP o kerma de referencia a la magnitud usada por el modelo.
- `table_normal` y `table_angular`, o alternativamente `table_total`: transmisión/atenuación de la mesa. No deben declararse ambas formas, porque producirían doble conteo.
- `bsf`: factor de retrodispersión.
- `f_skin`: conversión de kerma en aire a la magnitud de dosis en piel.

Para la etapa `skin_dose` deben existir `cf`, un modelo de mesa y `bsf`+`f_skin`, salvo que se utilice `surface_factors`.

## 6. Proveedor tabulado Benmakhlouf

La configuración de investigación usa:

```json
"dose_medium": "water",
"surface_factors": {
  "provider": "benmakhlouf_2011",
  "dataset_directory": "../datos/benmakhlouf_2011",
  "inherent_al_mm": 3.1,
  "target": "W",
  "field_shape_policy": "square_only",
  "out_of_domain": "error",
  "incident_kerma_basis": "free_in_air_without_backscatter",
  "quality_basis": "event_kvp_filtration_proxy_at_surface",
  "allow_water_phantom_approximation": true
}
```

En este modo, las tablas versionadas de `datos/benmakhlouf_2011` proporcionan `k_BS` y `k_f`. Se interpolan dentro del rango de sus nodos; si el kVp, la filtración, el campo o la geometría están fuera del dominio, la ejecución se detiene con un error controlado. El bloque `surface_factors` no debe combinarse con `bsf` o `f_skin` históricos.

## 7. `tube_output` y opciones adicionales

`tube_output` permite convertir mAs mediante un rendimiento documentado:

```json
"tube_output": {
  "reference_distance_mm": 1000,
  "yield_gy_per_mas": {
    "value": 0.0003,
    "provenance": "calibración",
    "identifier": "YIELD-v1"
  }
}
```

`rdsr_relative_standard_uncertainty` registra la incertidumbre asociada al RDSR. `allow_uniform_sweep` solo debe activarse cuando se haya documentado que un evento representa un barrido uniforme.

## 8. Qué no debe hacerse

- No cambiar el fabricante/modelo del JSON para ocultar una incompatibilidad del RDSR.
- No modificar el DICOM para rellenar campos ausentes.
- No usar a la vez `table_total` y `table_normal`/`table_angular`.
- No usar `surface_factors` junto con `bsf` o `f_skin`.
- No aceptar extrapolaciones de las tablas de Benmakhlouf.
- No tratar una configuración sintética como validación clínica.

Los resultados incluyen una auditoría por evento, el proveedor de factores, los nodos utilizados y el estado de interpolación. Esa información debe conservarse junto con `resultados.json` para que la ejecución sea reproducible.

## 9. Plantilla recomendada

La plantilla más completa disponible en el proyecto es:

`examples/equipo_xcat_benmakhlouf_investigacion.json`

Para una demostración sintética puede utilizarse:

`datos/equipo_sintetico.json`

La primera está orientada a XCAT + Benmakhlouf y la segunda a una superficie plana con factores sintéticos históricos.
