# Factores superficiales de Benmakhlouf: implementación y verificación

Fecha: 6 de octubre de 2026. Dataset `KS-ASF-201101-IR`, versión 1.0.0.
Estado: investigación; verificación numérica y del formalismo, **no validación clínica**.

## Fuente y análisis físico

La fuente primaria es el informe adjunto *Backscatter factors and mass energy-absorption coefficient ratios for surface dose determination in diagnostic radiology*, KS-ASF-201101-IR, Karolinska University Hospital, 7 de abril de 2011, de **Hamza Benmakhlouf, Annette Fransson y Pedro Andreo**. La portada de este PDF no incluye a Bouchard. Esta implementación no intercambia sus tablas con las de una versión de revista. El SHA-256 del PDF es `0d8ced5181445475d8540fd296d917d020f485ef75d61303003d425312d3f896`.

Las ecuaciones (1), (8) y (9) definen dos magnitudes adimensionales:

\[
k_{BS}=B_{air}=\frac{(K_{air})_{w,\mathrm{sup}}^{P+S}}{(K_{air})_{air}^{P}},\qquad
k_f=\frac{\int E\,\Phi_{\mathrm{sup}}^{P+S}(E)(\mu_{en}/\rho)_w(E)\,dE}
{\int E\,\Phi_{\mathrm{sup}}^{P+S}(E)(\mu_{en}/\rho)_{air}(E)\,dE}.
\]

Por tanto, ecuaciones (2) y (5):

\[
(K_w)_{w,\mathrm{sup}}=K_{air,\mathrm{inc}}\,k_{BS}\,k_f\simeq D_{w,\mathrm{sup}}.
\]

El cociente es **agua/aire para el espectro total en superficie**, no un coeficiente másico aislado, ni el cociente piel/aire, ni un cociente evaluado sólo con el espectro incidente. El producto transforma kerma en aire incidente libre en aire, en el mismo punto geométrico, en kerma superficial en agua. La aproximación a dosis supone depósito local de energía de electrones y pérdidas radiativas despreciables. Ambos factores dependen de calidad y campo. En el formalismo radioterápico, en cambio, `B_w` se combina con el cociente del espectro primario: no se mezclan ambas parejas.

Se emplean PENELOPE 2008/PenEasy, corte fotónico de 1 keV y electrones sin transporte explícito. La sección 3.2 describe agua de 15 cm de espesor, FSD de 100 cm, vóxel central de 1 × 1 × 0,01 cm y extensión lateral suficiente para dejar más de 5 cm de margen a la salida. Los campos son cuadrados de 5, 10, 20, 25 y 35 cm de lado en la superficie. Las calidades clínicas se generan con SpekCalc para tungsteno; mamografía utiliza Boone et al. El informe no proporciona un ángulo de ánodo numérico para cada fila clínica.

Hay inconsistencias que se conservan: los pies de tablas 1 y 3 indican 15 × 30 × 30 cm, incompatible con el mayor campo y el margen lateral descrito; la tabla 3 imprime `h > 1` en 24 calidades con Al=3,1 mm y 120–150 kV, pese a definirlo como HVL1/HVL2. Estos valores no se corrigen ni se usan para interpolar. Los resultados que los utilizan como metadato llevan `source_quality_flags`. El dominio de energías simulado en el texto es 4–150 keV (44 energías), pero la tabla 1 sólo publica 30 energías de 5–150 keV: no se reconstruyen las restantes.

## Datos extraídos y trazabilidad

Todos los archivos están en `src/skindose_verifier/data/benmakhlouf_2011/`.

| Archivo | Origen | Contenido y uso |
|---|---|---|
| `table1_monoenergetic.csv` | Tabla 1, p. 25 | 150 pares: 30 energías × 5 campos; consulta monoenergética independiente |
| `table2_comparison.csv` | Tabla 2, p. 26 | Seis espectros; BSF y cociente por MC completo y promedio, conservando los ratios impresos; comparación, no malla clínica |
| `table3_clinical.csv` | Tabla 3, pp. 27–32 | 715 pares: 143 calidades × 5 campos; proveedor para eventos |
| `table4_mammography.csv` | Tabla 4, p. 33 | 30 pares: seis calidades IEC × 5 campos; archivada, no activada en el proveedor intervencionista |
| `metadata.json` | Informe y análisis | Referencia, hash del PDF/CSV, versión, unidades, geometría, dominios y anomalías |
| `source_tables_layout.txt` | Páginas 25–33 | Extracción textual conservada para revisar las columnas originales |

Los decimales publicados se conservan como texto en CSV, incluidos ceros finales, HVL y `h`. Cada par tiene identificador de fila, tabla y página; en tablas 1/3 también fila de la página. Se preservan 1.814 valores de factores, contando las dos técnicas de tabla 2. La extracción reproducible con Poppler se contrastó mediante coordenadas de palabras con PyMuPDF: **2.577 celdas numéricas comparadas, ninguna discrepancia**. El control incluye todas las celdas numéricas de tablas 1/3 y todos los factores/ratios de tabla 2 y HVL/factores de tabla 4. Véase `BENMAKHLOUF_EXTRACTION_CHECK.json`.

## Selección e interpolación

`factor_tables.py` implementa un interpolador genérico de pares sobre mallas cartesianas completas. `surface_factors.py` carga datos y metadatos, comprueba hashes y selecciona la familia. El cálculo dosimétrico permanece en `engine.py`.

En tabla 3, cada familia de aluminio tiene una malla rectangular propia en kVp, Cu y lado del campo. Se fija el aluminio a una familia publicada, sin interpolar entre familias:

| Al inherente (mm) | Cu añadido tabulado (mm) |
|---|---|
| 2,5 | 0; 0,3; 0,9 |
| 3,1 | 0; 0,1; 0,2; 0,3; 0,6; 0,9 |
| 4,1 | 0; 0,1; 0,2 |
| 5,1 | 0 |

El voltaje abarca 50–150 kV en pasos de 10 kV. La selección usa kVp y Cu del evento y Al/ánodo declarados en configuración. HVL se conserva como descriptor; no se identifica un espectro por HVL sola, ni se reemplaza energía por kVp. La configuración exige reconocer que estas variables son un **proxy del espectro superficial**, no una caracterización espectral completa.

Se interpola linealmente en cada eje dentro de su intervalo. Según los ejes intermedios, resulta interpolación lineal, bilineal o trilineal; tabla 1 permite bilineal energía–lado. Los pesos son productos de fracciones lineales de distancia, positivos y de suma uno. Esto preserva los nodos, evita sobreoscilación de splines y respeta las dimensiones reales de las mallas. El eje de tamaño es **lado en cm, no área**. No se presume que la interpolación sea un nuevo cálculo Monte Carlo ni se le atribuye precisión clínica. Se interpolan ambos factores por separado y luego se multiplican; esto no equivale, en general, a interpolar su producto o a realizar una convolución espectral.

Fuera del dominio se lanza `InputError` con el UID del evento y el eje afectado: no se recorta, aproxima al vecino ni extrapola. La familia Al=5,1 sólo acepta Cu=0. Un campo rectangular se rechaza salvo opción explícita `equal_area`, que usa √área y registra una aproximación no validada por el informe.

## Integración y prevención de doble conteo

Antes de modificarlo, `engine._dose` multiplicaba kerma corregido por `bsf*f_skin`; `Factor` admitía constantes o tablas arbitrarias. La memoria existente, sección sobre `f_skin`, y `equipo_sintetico.json` muestran que el valor 1,05 era sintético y que `dose_medium` era una etiqueta, no una conversión material. **No puede demostrarse una equivalencia física general de ese factor con el de Benmakhlouf.** Se conserva íntegra la vía anterior y se añade un proveedor alternativo de pareja coherente. Sólo se permite con `dose_medium="water"`; se rechaza cualquier `bsf` o `f_skin` simultáneo, incluso en overrides por modo.

La ecuación aplicada por evento e y punto p es:

\[
D_w(p)=\sum_e M_e(p)K_{a,r,e}
\left(\frac{d_{r,e}}{r_e(p)}\right)^2
CF_e\,P_e(p)\,T_e(p)\,k_{BS,e}(p)\,k_{f,e}(p).
\]

`M` es la máscara geométrica de primera entrada visible, `CF` la calibración del indicador, `P` el perfil opcional y `T` la transmisión de mesa cuando se cruza, o uno en otro caso. El inverso del cuadrado se aplica una vez. `CF`, `P` y `T` deben referirse a kerma incidente y no incorporar retrodispersión ni conversión de medio; la configuración exige declarar esa base física, aunque no puede verificar cómo se midió una calibración externa.

El campo se escala a la profundidad axial de cada punto, en el plano perpendicular al haz, utilizando la referencia detector/isocentro existente. No se usa como tamaño de campo el área cutánea oblicua ni el área de una cara triangular. La distancia radial al foco se usa para kerma y la profundidad axial para colimación, sin confundir ambas.

Cada muestra angular registra `surface_factor_audit_per_sample`: consulta, familia, nodos originales, ecuación de pesos, estado tabulado/interpolado, hash, referencia, aproximaciones y anomalías. Los rangos `k_BS`/`k_f` son alias de auditoría de `bsf`/`f_skin`, **no multiplicaciones adicionales**. La incertidumbre de esta pareja queda no cuantificada: no se presenta el 0,3% del informe como incertidumbre total del mapa ni se suman como independientes factores correlacionados.

## Verificación independiente y tests

`validate_benmakhlouf.py` usa `Decimal` y una interpolación sucesiva propia sobre CSV para calcular expectativas; no llama al interpolador productivo para obtenerlas. Compara siete condiciones (nodos e interiores), y para cinco condiciones clínicas compara además el mapa con la ecuación de fuente puntual, geometría sintética y T=0,8.

| Condición (Al=3,1 mm para haces clínicos) | kBS | kf | Producto |
|---|---:|---:|---:|
| 50 kV, Cu=0, lado 5 cm | 1,204 | 1,016 | 1,223264 |
| 80 kV, Cu=0,1, lado 20 cm | 1,476 | 1,024 | 1,511424 |
| 150 kV, Cu=0,9, lado 35 cm | 1,592 | 1,067 | 1,698664 |
| 85 kV, Cu=0,15, lado 15 cm | 1,458625 | 1,02825 | 1,49983115625 |
| 73 kV, Cu=0,24, lado 12 cm | 1,44006 | 1,025584 | 1,47690249504 |
| Monoenergético 40 keV, lado 10 cm | 1,403 | 1,016 | 1,425448 |
| Monoenergético 42,5 keV, lado 15 cm | 1,4655 | 1,018 | 1,491879 |

Máxima diferencia absoluta de factores/producto: **2,22 × 10⁻¹⁶**. Máxima diferencia de dosis: **1,04 × 10⁻¹⁷ Gy**; redondeo de coma flotante, sin discrepancias de implementación identificadas. Resultados completos: `BENMAKHLOUF_VALIDATION.json`.

La tabla 2 del propio informe muestra diferencias entre promedio y MC completo de hasta 0,753% en BSF, 1,271% en kf y 0,753% en el producto, usando sus números redondeados. Se documentan como diferencias de método en la fuente, no errores del software; no se corrigen las tablas para eliminarlas. La expresión de promedio de la ecuación (12) no se ha utilizado para reconstruir espectros: se conserva la definición de la ecuación (9) y los valores publicados.

Se añadieron **33 casos de test** en `tests/test_surface_factors.py`: igualdad exacta de los 865 pares de tablas 1/3; interpolación manual lineal/bilineal/trilineal; extremos, NaN e infinito; dominios por familia; rechazo de doble conteo y medio incompatible; archivos propios/hash; acumulación por evento; correspondencia con dosis final y política rectangular. Suite completa ejecutada: **87 tests superados** (54 existentes + 33 nuevos), sin fallos.

## Uso y actualización

`examples/equipo_benmakhlouf_investigacion.json` es un ejemplo completo; mantiene geometría, CF y mesa **sintéticos**, pero activa datos publicados en agua. No es una calibración del equipo del hospital. Sobre un RDSR sintético generado con `demo`, ejecutar:

```bash
python -m skindose_verifier demo --out demo_sintetico --stage free_air
python -m skindose_verifier run demo_sintetico/rdsr_sintetico.dcm --config examples/equipo_benmakhlouf_investigacion.json --out salida_benmakhlouf --stage skin_dose
python -m pytest
python scripts/validate_benmakhlouf.py
```

La consulta directa admite, por ejemplo:

```python
from skindose_verifier.surface_factors import benmakhlouf_2011
r = benmakhlouf_2011().evaluate("clinical", inherent_al_mm=3.1,
                              kvp=80, added_cu_mm=0.1, field_side_cm=20)
print(r.k_bs, r.k_f, r.audit)
```

Para datos experimentales futuros usar `provider="tabulated_water"` y `dataset_directory`, con otro identificador, versión, procedencia, hashes y mallas completas según el esquema. El proveedor no altera los CSV publicados. Deben medir la misma pareja física; factores compuestos o de otro material requieren otro adaptador explícito. La propagación existente de escalas de `bsf`/`f_skin` no está habilitada para este proveedor y falla explícitamente al no existir esos factores configurados; no modifica los valores publicados.

Cambios existentes: `config.py`, `engine.py`, `pyproject.toml` y `README.md`. Archivos nuevos: `factor_tables.py`, `surface_factors.py`, los seis archivos del dataset, el ejemplo de equipo, tres scripts de extracción/comprobación/validación, `test_surface_factors.py`, este informe y los dos JSON de resultados.

## Límites físicos restantes

El respaldo directo del informe alcanza los factores centrales de una superficie plana de agua en sus calidades, campos y geometría. Su reutilización local sobre una malla es una aproximación declarada. No modela composición anatómica, hueso/pulmón, espesor real, bordes del campo, ausencia de equilibrio lateral, oblicuidad, ni retrodispersión espacial fuera de la máscara primaria. Tampoco calcula endurecimiento espectral por mesa/colchón: multiplicar T escalar no modifica la selección espectral. La reutilización a otra FSD se apoya en la dependencia débil descrita por los autores, pero no valida todas las geometrías de pacientes.

No se obtiene dosis en piel real ni dosis glandular por cambiar una etiqueta; el resultado es una estimación de **dosis superficial equivalente en agua**. Faltan caracterización espectral del equipo, transmisión y calibración medidas, incertidumbre de interpolación y geometría, y validación experimental sobre fantomas adecuados antes de cualquier interpretación clínica.

Como reparación de la importación inicial, se restauró también `data/pyskindose_adult_female.stl`: el remoto contenía 786.443 bytes y no producía caras válidas; el original conserva 1.337.884 bytes, 26.756 caras y SHA-256 `f74681f580aa4b5f281d7f6d49eb9888ffd25c25c9b129f355bf01d0b299cd91`, coincidente con `docs/ASSET_ORIGEN.md`. No se ha cambiado la anatomía ni la licencia.
