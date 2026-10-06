# Calculador y verificador de dosis en piel

**Versión 0.5.3.** Implementación **de investigación** del documento «Sistema independiente de dosimetría en cardiología intervencionista». Lee RDSR DICOM, reconstruye eventos y geometría del haz, calcula mapas de kerma/dosis sobre un plano, un cilindro elíptico **o una malla triangular 3D de piel**, acumula eventos y registra una PSD puntual y una media local de 1 cm² cuando la discretización lo permite. Incluye una rama de comprobación de Ka,r mediante rendimiento del tubo y pruebas sintéticas reproducibles.

**Estado físico:** el ejemplo incluido es completamente ficticio. Ninguna constante de su fichero `equipo_sintetico.json` es una calibración clínica. El programa imprime `research_unvalidated_no_clinical_use` en los resultados. Para comparar con casos reales hacen falta RDSR anonimizados, caracterización del equipo, dosimetría experimental, definición exacta de la magnitud y aprobación del protocolo de validación. El programa no reemplaza la medida de referencia ni el juicio del radiofísico.

## Ejecutar en Windows con Visual Studio Code

1. Descomprime el proyecto, abre la carpeta `skin_dose_project` desde **Archivo → Abrir carpeta** en VS Code y abre **Terminal → Nuevo terminal**.
2. Escribe estos comandos en PowerShell, **uno tras otro**:

   ```powershell
   python --version
   python -m venv .venv
   .\.venv\Scripts\python.exe -m pip install -e ".[test]"
   .\.venv\Scripts\python.exe -m pytest
   .\.venv\Scripts\python.exe -m skindose_verifier demo --out demo_output --stage skin_dose
   ```

3. Abre `demo_output/mapa_3d_interactivo.html` con el navegador para girar el mapa y mover sus planos de corte; también se generan `mapa_2d.png`, `mapa_3d.png` y `resultados.json`. `mapa_dosis.npz` guarda los valores acumulados **y los de cada evento**. La instalación de dependencias necesita Internet una vez; después el HTML funciona sin conexión.

Si aparece «`py` no se reconoce», utiliza `python` como en los comandos de arriba o haz doble clic en `EJECUTAR_DEMO.bat`: prueba ambos nombres automáticamente. Si **tampoco** se reconoce `python`, instala Python 3.10 o posterior desde [la página oficial](https://www.python.org/downloads/) y cierra **por completo** VS Code antes de volver a abrirlo para que detecte el comando nuevo. El instalador puede ofrecer incorporar los comandos de Python al `PATH`; permite esa opción si aparece.

Para ver el **fantoma humano triangular XCAT** del repositorio indicado, en vez del plano:

   ```powershell
   .\.venv\Scripts\python.exe -m skindose_verifier demo --model mesh --mesh-edge-mm 10 --out demo_mesh --stage skin_dose
   ```

Abre `demo_mesh/mapa_3d_interactivo.html` con doble clic: arrastra para girar, usa la rueda para acercar, selecciona dosis acumulada o un evento, y mueve los planos de corte X/Y/Z. Para explorar el ejemplo incluido **sin instalar Python**, ejecuta `VER_MAPA_XCAT.bat` o abre `results/research_xcat/mapa_3d_interactivo.html`. El STL original, la licencia y la procedencia están en `skindose_verifier/data/` y `ASSET_ORIGEN.md`.

También puedes abrir `EJECUTAR_DEMO.bat` con doble clic si ya tienes Python 3 instalado. No requiere Perl, MiKTeX ni compilar LaTeX: es un programa en Python independiente del PDF que describe el proyecto.

En Linux/macOS, crea un entorno con `python3 -m venv .venv` y utiliza `.venv/bin/python` en los comandos anteriores. Requiere Python ≥ 3.10.

## Comandos

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier inspect "archivo_anonimizado.dcm"
.\.venv\Scripts\python.exe -m skindose_verifier demo --out demo_elipse --model ellipse --stage skin_dose
.\.venv\Scripts\python.exe -m skindose_verifier demo --out demo_mesh --model mesh --stage skin_dose
.\.venv\Scripts\python.exe -m skindose_verifier run "archivo_anonimizado.dcm" --config "equipo_validado.json" --out "salida" --stage free_air
.\.venv\Scripts\python.exe -m skindose_verifier convergence "archivo_anonimizado.dcm" --config "equipo_validado.json" --stage skin_dose --resolutions 10,5,2.5
.\.venv\Scripts\python.exe -m skindose_verifier sensitivity-position "archivo_anonimizado.dcm" --config "equipo_validado.json" --offset-mm 5 --stage skin_dose --out "sensibilidad_posicion.json"
.\.venv\Scripts\python.exe -m skindose_verifier validate-measurements "mediciones.csv" --stage skin_dose
.\.venv\Scripts\python.exe -m skindose_verifier validate-map "salida/resultados.json" "mapa_medido.csv" --threshold 0.5 --max-distance-mm 5
.\.venv\Scripts\python.exe -m skindose_verifier compare-radimetrics "salida/resultados.json" --kar 0.1 --psd 0.15 --medium skin
```

`inspect` devuelve solo códigos DICOM, unidades, tipos y cantidades; no muestra identificadores personales. `run` necesita una configuración específica que concuerde exactamente con el fabricante/modelo del RDSR. Las tres etapas son `free_air` (kerma en aire libre a la superficie), `calibrated_air` (kerma en aire incidente tras CF y mesa) y `skin_dose` (dosis absorbida estimada en el medio declarado tras BSF y `f_skin`). **No compares valores de etapas o medios distintos.**

## Configuración de un equipo

El único fichero listo para ejecutar, `skindose_verifier/data/equipo_sintetico.json`, describe una simulación. Hay que elaborar otro JSON por equipo y versión con su fabricante/modelo, procedencia de calibraciones, dimensiones y convenciones geométricas. No copies los factores de ejemplo. Los parámetros clave son:

| Clave | Significado y unidad |
| --- | --- |
| `source_to_isocenter_mm`, `source_to_detector_mm`, `source_to_irp_mm` | Distancias de referencia; se usan como respaldo explícito si faltan en el evento y queda constancia. |
| `field_reference` | `detector` o `isocenter`; plano donde se reporta el ancho/alto del campo en el RDSR de ese equipo. |
| `primary_sign`, `secondary_sign`, `rotation_order` | Convención medida de ángulos y composición de giros. El ejemplo toma +y craneal, +z anterior y +x izquierda del paciente. |
| `table_coordinates`, `table_origin_mm`, `table_matrix` | Origen y transformación verificados de la posición de mesa RDSR a desplazamiento en coordenadas de paciente. |
| `table_plane_z_mm`, `table_half_width_mm`, `table_half_length_mm` | Aproximación geométrica de mesa y colchoneta; la transmisión solo se aplica si el rayo cruza su superficie acotada. |
| `patient` | `plane`, `ellipse` o `mesh`, dimensiones o archivo STL/OBJ/PLY, desplazamiento sobre mesa y resolución; deben contrastarse con el fantoma/paciente real. |
| `factors` | CF, transmisión `table_normal`×`table_angular` **o** `table_total`, BSF y `f_skin`; cada uno con procedencia e identificador. |
| `factors_by_mode` | Opcional: calibraciones específicas por código DICOM del modo de irradiación. Al activarlo se bloquean los modos sin entrada. Se mezcla cada entrada con los factores generales. |
| `dose_medium` | `skin` o `water`, según la magnitud experimental definida, solo para etapa `skin_dose`. |
| `tube_output` | Rendimiento medido en Gy/mAs y distancia de referencia, potencialmente tabulado por `kvp`/`filtration_mm_cu`. Solo constituye una estimación independiente **si** la medida y el mAs son fiables. |

Ejemplo de *forma* de una curva medida para reemplazar un factor constante (los números siguientes son **ficticios**):

```json
"cf": {
  "provenance": "medida local, dosímetro/certificado y fecha DOCUMENTADOS",
  "identifier": "CF-2026-v1",
  "axes": {"kvp": [60, 100], "filtration_mm_cu": [0.1, 0.3]},
  "grid": [[1.01, 1.03], [0.99, 1.02]],
  "relative_standard_uncertainty": 0.02
}
```

El orden de dimensiones de `grid` es el de `axes`; no se permite extrapolar. Para seleccionar CF por modo, por ejemplo, añade `"factors_by_mode": {"44491008": {"cf": { ... }}}` y proporciona entradas para todos los modos que vayas a calcular. Si el material del filtro no es Cu, o hay filtros combinados que requieren un modelo espectral, el parser bloquea el cálculo. Si un campo geométrico crítico está ausente y no está definido como respaldo en la configuración, también se bloquea.

### Malla 3D de piel

El ejemplo `--model mesh` utiliza `WholeBody_Body.stl` del fantoma [mesh50_XCAT](https://github.com/BenAuer2021/Mesh-based-Human-Phantom-for-Simulation) **Whole Body**, un modelo masculino antropomórfico compartido para simulación. Se incluye solo su contorno corporal externo de 12 148 triángulos; no se calculan dosis de órganos. La unidad mm, la matriz de ejes, la traslación aproximada a la mesa y la procedencia se documentan en `ASSET_ORIGEN.md`. El modelo MakeHuman anterior de PySkinDose sigue disponible con `--model mesh-pyskindose`.

Para usar otra malla de superficie cutánea cerrada (STL/OBJ/PLY), configura `"type": "mesh"`, `"mesh_path": "mi_piel.stl"` (ruta relativa al JSON), `"mesh_unit": "mm"` o `"cm"`/`"m"`, `"mesh_axes_matrix"` (matriz ortonormal 3×3), `"mesh_origin_mm"` (traslación tridimensional), `"origin_on_table_mm"`, `"thickness_mm"` y `"resolution_mm"` nominal. El programa verifica que la malla sea cerrada, conexa y tenga normales coherentes; comprueba la **primera cara intersectada** por cada rayo y suma la dosis por triángulo en el marco del paciente. Guarda vértices, caras, centros, áreas y dosis en `mapa_dosis.npz`.

El tamaño real de cada cara determina la resolución. `resolution_mm` sigue siendo nominal para mallas externas; ahora `mesh_max_edge_mm` subdivide realmente las aristas. `convergence` compara mallas refinadas con los tamaños indicados. Opcionalmente `mesh_refinement_roi_mm: [[xmin,ymin,zmin],[xmax,ymax,zmax]]` limita la región que solicita refinamiento, manteniendo sus vecinos conformantes. El límite por defecto es un millón de caras y el refinamiento nunca suaviza la anatomía. La media de 1 cm² puede ser `null` si el mallado local no define esa área con fiabilidad. No se debe extrapolar la forma de una persona a otra sin medidas anatómicas y validación.

### Visor 3D y comprobaciones espaciales

`mapa_3d_interactivo.html` es un archivo autónomo: se abre en un navegador con gráficos 3D activados y no descarga nada. Presenta la **misma escala de color** al alternar entre dosis acumulada y eventos. La superficie se puede girar y acercar; las vistas posterior, anterior, lateral y superior se activan con botones. Los planos X, Y y Z se desplazan con sus deslizadores y pueden dejar visible cualquiera de los dos lados. Su intersección recorta gráficamente cada triángulo; **la dosis original de la cara no cambia ni se interpola**. En los modelos plano y elíptico, muestra los puntos de muestreo sin inventar caras nuevas. La visualización del navegador puede utilizar la GPU, pero el **cálculo de dosis sigue en CPU** y no es en tiempo real durante un procedimiento.

En `resultados.json`, `surface_dose_area_approx` suma las áreas representadas por celdas cuya dosis supera 0, 10 %, 50 % o 90 % del máximo puntual. Son medidas descriptivas **aproximadas**, afectadas por la resolución y el tamaño desigual de los triángulos; no son límites clínicos ni histogramas de órganos. El array `event_dose_gy` de `mapa_dosis.npz` conserva cada evento en el mismo orden de `event_uid` y su suma reproduce exactamente `dose_gy`. Además, `global_psd_event_contributions` ordena los eventos según su contribución en el punto de PSD acumulada, lo que permite distinguir un campo dominante de una superposición de campos.

Para localizar el RDSR dentro de una carpeta DICOM exportada desde PACS o consola del equipo:

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier find-rdsr "datos\procedimiento_01"
```

El comando lee solo metadatos, lista los candidatos `X-Ray Radiation Dose SR`, resume los DICOM que no son RDSR, como imágenes `XA`, y devuelve hashes de UIDs en lugar de identificadores directos. Después se puede pasar la ruta candidata a `inspect` o `run`.

La salida `resultados.json` incluye también `factor_metadata`, con procedencia, identificador, representación, incertidumbre relativa y dominio tabulado de cada factor de corrección. Su objetivo es distinguir factores sintéticos de demostración frente a factores medidos, publicados o generados por simulación.

Para un pre-chequeo más completo de un RDSR real anonimizado:

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier compatibility-report "datos\procedimiento_01" --config "equipo_validado.json" --stage skin_dose
```

Este informe indica si el parser puede leer el RDSR, si faltan campos críticos por evento, si hay calidad de haz suficiente para la etapa dosimétrica solicitada y si fabricante/modelo coinciden con la configuración. No sustituye a una validación física.

`sensitivity-position` recalcula seis escenarios desplazando el paciente ±5 mm (o el valor elegido) en X, Y y Z respecto a la mesa. Presenta el rango de PSD y la variación de la posición del máximo. Solo mide la **sensibilidad bajo esos desplazamientos hipotéticos**; no proporciona un intervalo de confianza ni sustituye medidas geométricas del equipo. Los motivos bibliográficos y límites se detallan en `MEJORAS_Y_REFERENCIAS.md`.

## Medidas y comparación espacial

El CSV de medidas puntuales lleva cabecera `test_id,rdsr_path,config_path,measured_max_gy,standard_uncertainty_gy`; las rutas son relativas a la ubicación del CSV. **La medida debe corresponder a la magnitud y área de PSD del cálculo.** La salida de `validate-measurements` mantiene `validated: false` hasta que un profesional documente el protocolo, presupuesto de incertidumbre y criterios de aceptación.

El CSV de comparación espacial lleva `x_mm,y_mm,z_mm,measured_gy,standard_uncertainty_gy,quantity`, con al menos cuatro puntos en coordenadas locales del paciente. `quantity` debe ser literalmente la etiqueta de `resultados.json`, p. ej. `entrance_absorbed_dose_in_skin_Gy`. La herramienta compara solo posiciones medidas que coincidan con la malla dentro de la tolerancia, y devuelve RMSE, separación de picos y coincidencia de puntos sobre una isodosis. No registra automáticamente imágenes ni convierte la huella del dosímetro.

La comparación con Radimetrics es secundaria: la PSD debe estar definida en el mismo medio; la diferencia de Ka,r reportado **no es una verificación independiente** de su indicador. Conserva la rama de rendimiento de tubo para ello.

## Recorrido del documento de diseño

| Fases | Implementación y evidencia disponible |
| --- | --- |
| 0–2 Inventario, lector y normalización | Inventario de códigos/unidades, lectura de eventos y total acumulado, controles de ausencia, esquema por equipo y por modo. Probado solo con RDSR sintéticos; quedan pendientes RDSR reales y convención de fabricante. |
| 3–8 Coordenadas y geometría | Transformación de mesa configurada, dos ángulos, foco/detector, haz divergente, plano/elipse/malla corporal XCAT, primera intersección, distancia y área de campo. Tests sintéticos. Quedan pendientes fantoma físico y posicionamiento real. |
| 9–10 Dosimetría y coeficientes | Etapas dosimétricas separadas, interpolación acotada CF/mesa/BSF/`f_skin`, soporte por modo, estimación de incertidumbre **parcial**. Quedan pendientes curvas calibradas, dosímetros y medidas. |
| 11–13 Mapa, PSD y IRP | Suma por evento sobre la misma piel, máximos puntuales y vecindad de 1 cm², mapas por evento, superficie irradiada aproximada, salida 2D/3D e HTML interactivo, rendimiento de tubo condicionado, KAP/área etiquetado como consistencia. |
| 14–15 Validación e incertidumbre | Pruebas sintéticas, comparación de máximos medidos y puntos de mapa, convergencia de malla, análisis opcional de desplazamientos del paciente, lista de términos de incertidumbre pendientes. Las validaciones experimentales y clínica siguen pendientes. |
| 16–20 Software, trazabilidad y comparación | Paquete Python, comandos, tests, hashes de entrada/configuración, versiones de factores, comparación exploratoria con Radimetrics. No hay serie clínica suministrada. |

## Alcance geométrico y seguridad de interpretación

El campo se representa como rectángulo proyectado desde el foco, uniforme salvo perfil calibrado; el modelo elíptico es un cilindro **sin tapas**; el modelo de malla usa una geometría antropomórfica **sin órganos internos** y resolución heterogénea. La mesa es un plano acotado, y el paciente se considera rígido. La angulación durante un evento se bloquea salvo si se declara explícitamente `allow_uniform_sweep=true` y `sweep_samples>=2`; la dosis uniforme durante el barrido es entonces una **hipótesis**. No se infiere movimiento del paciente, perfil temporal, efecto talón ni transporte Monte Carlo. Solo se permite representar un perfil transversal mediante `beam_profile` calibrado y explícito; no se modelan físicamente filtros de cuña. Un presupuesto parcial de incertidumbre no equivale a un intervalo clínico.

Se guardan el hash SHA-256 de entrada/configuración, la versión de software/configuración, las versiones de factores, estados de calidad, resultados por evento y advertencias de acumulados. El RDSR original no se copia a la salida de `run`.

Para repetir la verificación del paquete: `python -m pytest` y `python -m compileall -q skindose_verifier` dentro del entorno virtual. Los ensayos son de consistencia computacional; **no son validación física**.


## Fase de I+D v0.5: uso y documentación

El ejemplo incluido usa el mismo STL XCAT original con refinamiento global a aristas de 10 mm: 177 926 caras. El visor permite fuente/isocentro/haz, arco esquemático, elección de evento y posición de barrido, inspección por clic de centroides irradiados, coordenadas, normal, incidencia, distancia y contribuciones. Los colores suavizados son opcionales y solo gráficos. Usa **Acercar a la zona irradiada** para ver el detalle; los cortes no recalculan dosis. El inspector corresponde al centroide de la cara, no a cualquier punto continuo pulsado.

`surface_metrics` contiene Dmáx y Dmedia/D2%/D5%/D10% ponderados por área, tanto sobre toda la superficie modelada como sobre la irradiada. Esta distinción evita percentiles ambiguos en un cuerpo con mucha piel a dosis cero. `histograma_dosis_area.csv` conserva A(D≥umbral), incluido el área completa a umbral cero; la gráfica comienza en el área irradiada (0+). No hay DVH volumétrico porque no hay espesor de piel. Las regiones ≥90% se agrupan por aristas compartidas; no constituyen un umbral de toxicidad.

Para el estudio completo, ejecuta **EJECUTAR_INVESTIGACION.bat** después de instalar el entorno con `EJECUTAR_DEMO.bat`. Ejecuta pruebas, demo refinada, convergencia, benchmark y propagación con 64 y 256 muestras. Puede tardar varios minutos y consumir bastante más memoria con 709 122 caras a 5 mm. No requiere Perl ni MiKTeX.

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier mesh-quality --config results/research_xcat/equipo_sintetico.json --out calidad.json
.\.venv\Scripts\python.exe -m skindose_verifier convergence results/research_xcat/rdsr_sintetico.dcm --config results/research_xcat/equipo_sintetico.json --resolutions 40,20,10,5 --out convergencia.json
.\.venv\Scripts\python.exe -m skindose_verifier uncertainty results/research_xcat/rdsr_sintetico.dcm --config results/research_xcat/equipo_sintetico.json --spec examples/incertidumbre_sintetica.json --out incertidumbre
.\.venv\Scripts\python.exe -m skindose_verifier benchmark results/research_xcat/rdsr_sintetico.dcm --config results/research_xcat/equipo_sintetico.json --out benchmark.json
```

Las distribuciones en `examples/incertidumbre_sintetica.json` son hipotéticas. Se declaran procedencia, semilla, número de muestras y distribución (normal/uniforme; lognormal para escalas positivas). `standard_deviation` es desviación típica absoluta de cada parámetro, también en lognormal, no el parámetro sigma de su logaritmo. Las escalas tienen media 1 por defecto; las traslaciones/ángulos, media 0. Los sesgos son comunes a todos los eventos de cada muestra. El archivo NPZ guarda media, desviación y frecuencia de exceder un umbral en cada centroide. No se calcula un presupuesto completo ni una probabilidad clínica de lesión.

Se documentan versiones de Python/dependencias y hashes del código en los resultados. `requirements-research-lock.txt` captura la ejecución Linux/Python 3.12; `pyproject.toml` conserva los rangos compatibles generales. No se ha ejecutado esta versión en Windows real. Para reproducir las versiones probadas, usar Python 3.12 y `python -m pip install -r requirements-research-lock.txt` antes de instalar el paquete; otras plataformas pueden necesitar resolver versiones compatibles.

Documentación científica: [auditoría](AUDITORIA_TECNICA.md), [revisión y referencias](RESEARCH_REVIEW.md), [resultados de verificación](RESULTADOS_DE_VERIFICACION.md), [protocolo experimental](VALIDACION_EXPERIMENTAL.md), [trazabilidad](CHANGELOG_RESEARCH.md) y [roadmap](ROADMAP_CALCULADOR_DOSIS_PIEL.md).
