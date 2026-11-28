# Sistemas de referencia: del TC al isocentro

DermaDose coloca la malla del paciente siempre con la misma cadena de transformaciones rígidas: TC (C) → tablero (T) → isocentro (I). Cada eslabón puede proceder de una medida o de una estimación declarada. El cálculo de dosis no distingue entre ambas; los resultados indican qué eslabones se han estimado.

## Marcos

| Marco | Origen | Ejes positivos |
|---|---|---|
| **C · TC** | Cero de las coordenadas físicas DICOM del TC | x izquierda del paciente, y posterior, z cabeza (LPS) |
| **T · tablero** | Centro del borde de cabecera, sobre la cara superior del tablero | x izquierda del tablero, y abajo, z hacia la cabecera |
| **I · angiógrafo** | Isocentro del arco | Paralelos a T con el tablero nivelado y sin giro |

Una malla exportada en RAS se convierte a LPS mediante (x, y, z) → (−x, −y, z). Internamente, el motor usa los ejes (izquierda, cabecera, arriba) para los ángulos del arco. `ROOM_FROM_I` es la permutación fija entre I y esos ejes, sin traslación ni datos de sala.

## Ecuaciones

$$p_T = R_{T\leftarrow C}\,p_C + t_{T\leftarrow C},\qquad t_{T\leftarrow C} = a_T - R_{T\leftarrow C}\,a_C,\qquad a_T = (0,\,-e,\,-d)$$

$$t_k = t_0 + S\,(m_k - m_0),\qquad p_I(k) = R_{I\leftarrow T}\,p_T + t_k,\qquad R_{I\leftarrow T} = \mathrm{Id}_3$$

$$M_{I\leftarrow C}(k) = M_{I\leftarrow T}(k)\,M_{T\leftarrow C}$$

- `a_C`, `a_T`: el mismo anclaje en el TC y en el tablero. `e` es el espesor de la colchoneta bajo carga y `d` la distancia desde la cabecera.
- `m_k`: lecturas longitudinal, lateral y altura del evento k (DCM 113751, 113752 y 113753), en mm.
- `S`: matriz 3 × 3. Su columna j es el desplazamiento en I por mm de cambio de la lectura j.
- `t_0`: posición de O_T en I cuando las lecturas valen `m_0`.

El ejemplo numérico del documento de referencia, (0, 90, 0) y (45, 60, 100) mm, es una prueba automática: `tests/test_reference_chain.py`.

## Cada eslabón y su evidencia

| Eslabón | Medido o explícito | Estimación, si falta | Dónde queda |
|---|---|---|---|
| `M_T←C` | `registro_tc.json` de DermaDose Surface (anclaje de hombros), anclaje y rotación manuales o tres o más marcas | `R` según la postura codificada del RDSR. Punto diana del TC sobre la línea media del tablero, con la espalda (la superficie de apoyo) sobre la colchoneta, a una distancia `d` de la cabecera | `patient.ct_registration` (`status`, `estimate`) |
| `t_0`, `m_0` | Marcador con coordenadas en T y en I (medidas o dos vistas calibradas), con las lecturas de esa posición | Técnica isocéntrica: el punto diana estaba en el isocentro durante el evento de referencia (por defecto, el primero con kerma) | `geometry.table_calibration.reference` (`origin_isocenter_mm`, `status`, `estimate`) |
| `S` | `motion_matrix_i_per_mm` obtenida con movimientos de un solo eje | Ninguna. Una columna `null` solo se admite si esa lectura no cambia entre eventos | `table_calibration.motion_matrix_status` |
| `R_I←T` | Tablero nivelado y sin giro. Se comprueba con DCM 113754–113756 y con las posiciones finales si el RDSR las incluye | — | — |

El punto diana se elige en este orden:

1. El anclaje de DermaDose Surface que coincide con la región diana del RDSR.
2. El anclaje torácico de Slicer (punto medio pulmonar).
3. Un punto al nivel cardiaco, 180 mm caudal al anclaje posterior de hombros y a media profundidad.
4. El centro geométrico de la malla, marcado como tal.

El lanzador (`python -m skindose_verifier.example_runner CASO`) aplica esta tabla. Solo estima si el perfil declara `azurion_profile.allow_research_assumptions` o si se usa `--permitir-estimaciones`. `configuracion_efectiva.json` guarda los valores usados, de modo que el mismo cálculo se repite sin volver a estimar. En `resultados.json`, `placement_chain` resume la cadena y `ct_registration` contiene la evidencia completa.

Opciones del lanzador relacionadas: `--espesor-colchoneta-mm`, `--distancia-hombros-mm`, `--evento-referencia first|max_kerma|N`, `--registro-tc` y `--mesa`.

## Lecturas de mesa del Azurion

La declaración de conformidad DICOM del Azurion lista las tres lecturas en mm, sin convenio de signo, origen ni eje. La configuración de Radimetrics de la sala (hoja «Calibración Radimetrics», D57–D61) da los ejes. Los signos se fijaron con RDSR externos del Azurion en los que la mesa se desplazó una distancia conocida; el detalle está en la memoria.

| Lectura | Movimiento físico | Columna de `S` en el perfil | Evidencia |
|---|---|---|---|
| Longitudinal (113751) | Izquierda | (1, 0, 0) | Radimetrics D59 y desplazamientos de mesa en haces angulados a izquierda y derecha. |
| Lateral (113752) | Cabecera | (0, 0, 1) | Desplazamientos de mesa en haces craneales: la lectura baja cuando el centro del campo se mueve hacia los pies. |
| Altura (113753) | Arriba | (0, −1, 0) | Radimetrics D61 y una subida de mesa de 11 cm que lee 111 mm más. DICOM define el sentido contrario (hacia el suelo). |

En Philips, la lectura «longitudinal» mueve izquierda–derecha y la «lateral» cabeza–pies, al contrario que sus nombres. La geometría heredada anterior aplicaba el convenio contrario y se ha retirado de los perfiles Azurion.

Las ganancias coinciden con 1 mm/mm dentro de la tolerancia de centrado de la cámara (±1 cm). Una columna `null` sigue siendo posible en otros perfiles: un caso en el que esa lectura cambie se detiene con un mensaje. Para medir los signos y ganancias de otra mesa:

1. Haz las tomas 4–6 de «Calibración geométrica»: desde la frontal, mueve solo un eje cada vez y anota hacia dónde va el tablero visto desde el paciente supino (y el recorrido, si lo mides con una marca de sala).
2. Ejecuta:

```powershell
python -m skindose_verifier.table_calibration serie.dcm --perfil config/mesa_azurion_hemodinamica_1.json `
    --movimiento 4:izquierda:100 --movimiento 5:arriba --movimiento 6:cabecera:150 --out mesa_medida.json
```

`EVENTO` es el número del evento en el RDSR; se compara con el anterior, o con `BASE>EVENTO`. Sin recorrido medido se supone ganancia unitaria y queda registrado. El perfil nuevo no sobrescribe el anterior.

## Referencia de sala medida

Para sustituir la estimación de `t_0` por una medida, con la mesa nivelada y quieta:

1. Fija un marcador y mide su posición en T. Con la distancia cabecera–marcador `d`, la altura sobre el tablero `h` y el lateral con signo `l`, `marker_table_mm = (l, −h, −d)`.
2. Guarda las lecturas del RDSR en esa posición en `readings_long_lat_height_mm`.
3. Localiza el marcador en I con `marker_isocenter_mm` o con `calibrated_views`: dos o más vistas con su matriz de proyección I → píxel y el píxel observado.

Con estos datos, `table_placement` pasa de `estimated_reference` a `marker_reference`.

## Límites

La cadena es rígida: no representa la respiración, la deformación por el apoyo ni los cambios de postura entre el TC y el procedimiento. Una referencia estimada sitúa el mapa de forma plausible para un caso isocéntrico, pero no es una calibración. La posición del máximo debe interpretarse como condicionada a las hipótesis que figuran en `placement_chain` y en `ct_registration.room_reference_estimate`.
