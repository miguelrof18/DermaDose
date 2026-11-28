# Cambios

## 0.8.1 · noviembre de 2026

### Corrección angular de la transmisión de mesa

- Perfiles `azurion_benmakhlouf_spekpy`, `ejemplo/equipo_azurion`, `benmakhlouf_tabulado` y `equipo_benmakhlouf_investigacion`: `table_angular` deja de valer 1 y pasa a la tabla de lámina `SLAB-T0.80`, T(θ)/T(0) = 0,8^(1/cos θ − 1), con θ el ángulo de cada rayo respecto a la normal del tablero (0,97 a 30°, 0,91 a 45°; constante desde 80°). Es un cambio de datos: el código no cambia. Los perfiles sintéticos de autoconsistencia conservan 1.
- Con T₀ = 0,8 fijo, en los casos oblicuos RFM 2024 la diferencia crecía +2,5 % por cada 10° (+6,3 % a 45°). Con la tabla de lámina la tendencia desaparece (−0,3 % por cada 10°) y los casos a incidencia normal no cambian ([detalle](docs/verificacion_rfm2024.md#transmisión-de-mesa-en-los-perfiles)).
- Informe de verificación y cifra de 0.8.0 actualizados al código vigente (área del Canon referida al detector): media +3,0 %, desviación típica 3,3 %, cuatro casos fuera de ±6,5 % y ninguno fuera de ±13 %.

## 0.8.0 · noviembre de 2026

### Verificación con el conjunto RFM 2024

- Verificación con la metodología de Gañán Mora et al. (Rev Fis Med 2024;25(2):25–34) sobre sus 27 RDSR ([informe y modelo](docs/verificacion_rfm2024.md)). Todos los casos quedan dentro de ±13 %. La cifra publicada en 0.8.0 (media +2,9 %, tres casos fuera de ±6,5 %) correspondía al área del Canon referida a 1000 mm; con el código vigente, ver 0.8.1.
- `surface_factors.aluminium_interpolation = linear_between_families`: interpolación lineal auditada entre familias de Al adyacentes, sin extrapolación. Antes, un Al total no tabulado detenía el cálculo (caso 6 del Canon, 1,8 mm Al).
- Perfiles Azurion: la columna lateral de `S` pasa de `null` a +1 (hacia la cabecera). Lo establecen los desplazamientos de mesa de los casos 16–17. Los casos 13–15, 18–19 y 24–27 confirman los signos longitudinal y de altura con ganancia unitaria (±1 cm).

## 0.7.0 · noviembre de 2026

### Colocación del paciente: una sola cadena TC → tablero → isocentro

- La malla se coloca siempre con `p_I(k) = R_T←C p_C + t_T←C + t_0 + S (m_k − m_0)` ([Sistemas de referencia](docs/sistemas_de_referencia.md)). Desaparecen la ruta provisional heredada (`provisional_ct_to_table`, `table_room_offset_mm`) y las colocaciones específicas del XCAT y del TC sin registro (`ct_placement.py`).
- Cada eslabón es medido o estimado, con su hipótesis guardada en `configuracion_efectiva.json` y en `placement_chain`:
  - `M_T←C` estimado: postura codificada del RDSR, punto diana del TC, espalda sobre la colchoneta, distancia a la cabecera.
  - `t_0, m_0` estimados: punto diana en el isocentro en el evento de referencia (`--evento-referencia`).
- `S` admite lecturas sin caracterizar (`null`). Si esa lectura cambia, el cálculo se detiene en lugar de suponer un sentido. Nueva herramienta `table_calibration` que obtiene `S` de movimientos de un solo eje (tomas 4–6 de «Calibración geométrica»).
- Perfiles Azurion: se retira la geometría heredada (`table_matrix`, origen y altura de mesa), que asignaba las lecturas longitudinal/lateral al revés que la configuración Radimetrics de la sala.
- DermaDose Surface registra `anchor_definition` en `registro_tc.json`.

## 0.6.0 · noviembre de 2026

La física del modelo no cambia: mismas ecuaciones, convenciones angulares, tablas y factores. En los casos de regresión comparados, los mapas, las métricas y la media de 1 cm² son idénticos a los de 0.5.3.

### Correcciones

- Tablas espectrales derivadas (`derived_mrc200_nominal`, `derived_azurion_spekpy`, `derived_azurion_qc_conditional`): Git convertía sus finales de línea y el SHA-256 dejaba de coincidir fuera de Windows. Ahora se conservan byte a byte (`-text`).
- `mesh_convergence` llamaba a un método inexistente.
- La propagación de incertidumbre descartaba el desplazamiento del paciente con registro TC–tablero provisional.
- El panel de factores del visor 3D quedaba siempre vacío.
- La sustitución por factores históricos resolvía rutas desde la carpeta de trabajo.
- El presupuesto parcial de incertidumbre nombraba `bsf`/`f_skin` con el proveedor tabulado.
- La colocación automática leía la postura del texto del RDSR en lugar de sus códigos.

### Rendimiento y memoria

- Primera entrada por proyección central vectorizada (`numpy_central_projection_cpu`), con Trimesh/VTK como respaldo (`DERMADOSE_RAY_BACKEND=trimesh`).
- Interpolación tabulada vectorizada, idéntica bit a bit.
- Mapas por evento dispersos en memoria (`EventMaps`) y en el visor HTML.
- Caché de geometrías limitada por memoria (256 MB) en lugar de dos entradas fijas.
- La búsqueda de la media de 1 cm² termina en cuanto ningún candidato puede mejorarla.

### Salidas y uso

- `resultados.json` añade `max_air_kerma_gy` (etapas `free_air` y `calibrated_air`) o `max_modeled_surface_dose_gy` (`skin_dose`), `psd_point_definition` y `dose_medium`, según la propuesta de la memoria (sección 2.4.1). `psd_point_gy` se conserva.
- Tiempos por fase en `performance`: preparación, proyección, dosimetría y métricas. El resultado de `run` añade la exportación.
- Progreso por evento con su duración. Opción `--quiet` en `run` y `demo`.
- Avisos para claves desconocidas en el JSON de equipo; los nombres que empiezan por `_` se consideran comentarios.
- Unidades UCUM adicionales en el lector RDSR.

### Pruebas

- Nuevas pruebas del núcleo de visibilidad (contra fuerza bruta), de integridad de datasets y del contrato de salida.
