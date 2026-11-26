# Registro del TC en la mesa del Azurion de Hemodinámica 1

DermaDose puede usar las dimensiones del tablero en cualquier caso de esta sala y, cuando se dispone de una referencia fija y un registro del paciente, situar cada punto de su TC respecto al isocentro en cada evento. La posición del máximo de dosis se devuelve en las coordenadas físicas del TC.

## Datos disponibles de la sala

| Dato | Valor del Excel | Uso en el programa |
|---|---|---|
| Tablero, largo × ancho × espesor | 3160 × 500 × 30 mm | Volumen rectangular que intercepta los rayos entre el foco y la piel. |
| Dirección positiva del detector | Primaria a la izquierda, secundaria hacia la cabeza | Orienta la fuente y el detector con los ángulos del evento. |
| Offsets del asistente Radimetrics | Lateral 508 mm, longitudinal −70,5 mm, altura −732 mm | Se conservan con las lecturas de regla, imagen y haz que figuran en el Excel. El asistente de Philips emplea también una posición del haz. Para llevar estos offsets al origen de DermaDose falta conocer el mismo punto físico en ambos marcos. |
| Traslación positiva de la mesa | Longitudinal hacia la izquierda, altura hacia arriba, lateral sobre el eje cabeza–pies | Los dos sentidos explícitos se recogen en el perfil. La celda D58 dice «Cranial/Caudal» y deja pendiente el signo de esa columna. |
| Hombros, colchoneta y altura al suelo | 400, 65 y 895 mm, respectivamente | Se conservan los rótulos originales. El registro de cada paciente requiere su distancia cabecera–hombros y el espesor de la colchoneta bajo su carga, si se usa el anclaje anatómico. |

La hoja «Calibración geométrica» contiene la plantilla de las tomas, sin valores del marcador ni lecturas de referencia. El perfil distribuido permite usar la geometría medida con la colocación previa del caso. La colocación absoluta se activa al completar las medidas de sala y proporcionar el registro del TC. La transmisión dosimétrica se sigue tomando del perfil de factores seleccionado. [Perfil de Hemodinámica 1](../config/mesa_azurion_hemodinamica_1.json)

Hay un signo que merece comprobarse en la sala. En D61, «Up» describe la configuración de Radimetrics. DICOM define el aumento de 113753 hacia abajo. Para resolver cuál es la relación efectiva en este equipo, compara la lectura del RDSR antes y después de subir físicamente el tablero. [Definición de las lecturas de mesa en DICOM PS3.16](https://dicom.nema.org/medical/dicom/current/output/chtml/part16/chapter_D.html#DCM_113753)

## Qué se mide una vez y qué corresponde a cada paciente

1. **Una posición de referencia de la sala.** Con el tablero nivelado y sin giro, elige un marcador fijo. Anota sus coordenadas en el tablero (`marker_table_mm`) y las tres lecturas del RDSR en esa misma posición (`readings_long_lat_height_mm`). Introduce las coordenadas 3D del mismo marcador respecto al isocentro en `marker_isocenter_mm`. Como alternativa, `calibrated_views` acepta dos o más vistas, con la mesa inmóvil. Cada vista lleva su matriz de proyección de I a píxeles `projection_px_from_i_mm` (3 × 4, profundidad positiva delante del foco) y el píxel observado `pixel_uv` (columna, fila). La matriz incluye la geometría y el origen real en la imagen. DermaDose triangula la marca y registra el error cuadrático medio de reproyección en píxeles.
2. **Movimiento de la sala.** Completa `motion_matrix_i_per_mm` con los signos y escalas comprobados. También puedes introducir tres entradas en `axis_trials`. Cada una contiene `delta_long_lat_height_mm` y `delta_isocenter_mm`, medidos tras mover un solo eje. DermaDose divide el desplazamiento físico por el cambio de lectura para obtener la columna correspondiente. Para comprobarlo, desplaza la mesa a otra posición y compara el marcador predicho con el observado en escopia y cine.
3. **Colocación del paciente.** Selecciona el mismo punto anatómico en el TC y en el tablero. Si usas el plano de hombros, registra el desplazamiento lateral, la distancia desde la cabecera y la altura sobre la cara superior. Aporta también la orientación del paciente. Un único punto fija la traslación una vez conocida la rotación. Si la orientación es desconocida, el programa admite tres o más pares de marcas 3D no colineales.

## Cómo se calcula la posición

| Marco | Origen | Sentido positivo de los ejes x, y, z |
|---|---|---|
| C, TC físico en mm (LPS) | Origen físico del TC | Izquierda, posterior, cabeza. Una malla exportada en RAS se convierte mediante $(x,y,z)\mapsto(-x,-y,z)$. |
| T, tablero | Centro del borde de cabecera, en la cara superior | Izquierda, abajo, cabeza. El tablero se extiende hacia los pies. |
| I, angiógrafo | Isocentro | Izquierda, abajo, cabeza en el montaje de referencia. T e I tienen ejes paralelos con el tablero nivelado y sin giro. |

Para un vértice de coordenadas $p_C$, el anclaje del paciente tiene coordenadas $a_C$ y $a_T$ en los respectivos marcos. La matriz $R$ gira los ejes de C a T. Si la línea media queda centrada y el plano de hombros está a distancia $d$ de la cabecera sobre una colchoneta de espesor comprimido $e$, entonces $a_T=(0,-e,-d)$. El desplazamiento del origen del TC en T es $a_T-Ra_C$.

En la referencia, $q_T$ y $q_{I,0}$ son las coordenadas del marcador en T e I. Para el evento $k$, $m_k$ contiene las lecturas longitudinal, lateral y altura en mm (DCM 113751, 113752, 113753). $m_0$ es la misma terna en la referencia fija. Las columnas de $S$ convierten el cambio de cada lectura en un vector de desplazamiento físico en I. [Definiciones DICOM PS3.16](https://dicom.nema.org/medical/dicom/current/output/chtml/part16/chapter_D.html)

$$
p_T=R(p_C-a_C)+a_T,\qquad
t_k=q_{I,0}-q_T+S(m_k-m_0),\qquad
p_{I,k}=p_T+t_k .
$$

$$
M_{T\leftarrow C}=
\begin{pmatrix}R&a_T-Ra_C\\0&1\end{pmatrix},\qquad
M_{I\leftarrow T,k}=
\begin{pmatrix}\mathrm{Id}_3&t_k\\0&1\end{pmatrix},\qquad
M_{I\leftarrow C,k}=M_{I\leftarrow T,k}M_{T\leftarrow C}.
$$

Las matrices actúan sobre vectores columna homogéneos. Un punto tiene cuarta componente 1 y una dirección tiene cuarta componente 0. El visor usa la inversa de la misma transformación para dibujar foco, detector y campo sobre la malla.

## Cómo ejecutarlo

Para usar las dimensiones y los sentidos del detector con la colocación habitual, desde la raíz del repositorio:

```powershell
python -m skindose_verifier.example_runner ejemplos/caso --mesa config/mesa_azurion_hemodinamica_1.json
```

Sustituye `ejemplos/caso` por la carpeta del caso, que contiene el JSON de equipo, el RDSR, el STL y su JSON de Slicer. Para activar el registro completo, guarda una copia del perfil con las medidas de referencia y proporciona el registro específico del TC:

```powershell
python -m skindose_verifier.example_runner ejemplos/caso --mesa sala/mesa_medida.json --registro-tc ejemplos/caso/registro_tc.json
```

El lanzador de Windows también detecta `mesa_azurion_hemodinamica_1.json` y `registro_tc.json` junto al caso. En un JSON de equipo se pueden definir `geometry.table_calibration` y `patient.ct_registration` como objetos o rutas relativas a ese JSON. El STL y su archivo homónimo deben documentar el sistema RAS o LPS, los milímetros, el nombre del fichero y su SHA-256.

Si toda la referencia de sala está ausente, un registro TC–tablero válido se conserva y se aplica también en el modo provisional. La geometría heredada (origen, matriz y altura de mesa) sitúa el tablero respecto al angiógrafo. Con el perfil provisional `reference_table_from_first_event`, el lanzador puede estimar una traslación horizontal de todo el tablero a partir de un punto diana anatómico del TC y de un haz suficientemente transversal del RDSR. Si no hay una región reconocida, se puede introducir un punto RAS conocido del procedimiento o dejar la referencia sin ajustar. Esta estimación no modifica `M_T_from_C`, no mide la altura ni los signos de movimiento y no es validación independiente. La traslación queda en `geometry.table_room_offset_mm` y su procedencia en `geometry.table_room_reference_estimate`. El diagnóstico conservador de intersección con la caja de la malla se exporta a `diagnostico_geometria.json` antes de calcular dosis; si ningún campo puede cruzarla, el cálculo se detiene inmediatamente. El ajuste anatómico no se modifica. Los resultados identifican este modo como `ct_table_explicit_room_provisional` y guardan `M_T_from_C`. Una referencia de sala parcialmente cumplimentada o un registro de paciente inválido siguen produciendo un error. No es una calibración clínica.

Los UID de marco distintos (por ejemplo, tras anonimización) requieren verificar todos los IOP e IPP contra la matriz IJK→LPS, continuidad y ausencia de duplicados, además de PixelSpacing y dimensiones cuando se exportan. Los metadatos antiguos se aceptan si conservan la geometría suficiente; las comprobaciones no disponibles se documentan. La coherencia de la rejilla no acredita la identidad del paciente ni el origen de los UID. Se exige vincular el registro a la malla mediante SHA-256. La postura del intervencionismo se contrasta con el RDSR, sin volver a aplicar PatientPosition del TC. Si difiere del HFS manual, el lanzador permite confirmar una rotación nominal conservando las dos marcas; las posiciones variables se rechazan. La calibración de sala sigue siendo independiente.

Este es un ejemplo **ficticio** de `registro_tc.json` con anclaje y rotación. Sustituye el identificador y todos los datos geométricos por los del caso. La identidad representa aquí una orientación alineada.

```json
{
  "mesh_sha256": "SHA-256 del STL de este caso",
  "anchor_ct_lps_mm": [20, 110, -140],
  "anchor_table_mm": [0, -35, -280],
  "rotation_table_from_ct": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
}
```

También se puede identificar el TC mediante su `frame_of_reference_uid` si coincide con el de los metadatos del STL. Para registrar con marcas, sustituye anclajes y rotación por `fiducials_ct_lps_mm` y `fiducials_table_mm`. Introduce los puntos en el mismo orden. El informe devuelve el RMS del ajuste en mm.

## Qué comprobar en el resultado

En `resultados.json`, `table_placement` distingue la referencia fija del marcador (`marker_reference`) de la colocación previa (`legacy_placement`). En el modo registrado, `ct_registration` identifica la malla y el ajuste utilizado. Cada muestra geométrica de cada evento contiene `M_T_from_C`, `M_I_from_T`, `M_I_from_C`, su inversa `M_C_from_I` y las tres lecturas originales normalizadas a mm. Tanto el pico como los puntos y las normales de `mapa_dosis.npz` se expresan en TC LPS.

Este perfil corresponde al arco único o al plano A. El registro fijo requiere las tres lecturas de mesa caracterizadas. Si el RDSR sólo contiene posiciones relativas al isocentro de otra familia DICOM, hace falta caracterizar esa correspondencia. Si el RDSR declara inclinación, giro o posiciones inicial y final diferentes durante una exposición, DermaDose señala el evento. Cuando faltan esos campos opcionales, el cálculo presupone que el tablero permaneció nivelado y quieto durante ese evento.
