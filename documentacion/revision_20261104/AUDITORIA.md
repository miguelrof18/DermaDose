# Revisión física y de software de DermaDose — 4 de noviembre de 2026

## Alcance y estado comprobado

Base: `main` c59ed266b34c3553b82216c12f977034919bd0c4. Revisión de código, instrucciones, memoria canónica, configuraciones, pruebas, tablas y propuestas abiertas. No se dispone de una sesión gráfica de Slicer ni de nuevas medidas del hospital. Esta revisión no certifica la ausencia de otros errores.

| Versión | Estado observado | Funcionalidad y límites |
|---|---|---|
| main c59ed26 | Integrada antes de esta revisión | Motor, RDSR, Benmakhlouf, SpekPy y módulo Slicer síncrono. |
| PR #9, feature/slicer-ct-surface | Borrador antiguo | Propuesta inicial de Slicer; parte de su funcionalidad ya existe en main. No debe confundirse con la versión asíncrona. |
| PR #10, codex/philips-azurion-profile | Borrador con conflictos, 14 commits detrás de main | Parte espectral ya está en main. Su descripción conserva políticas antiguas de incompatibilidad/fallback. No integrarla a ciegas ni usarla para sustituir main. |
| PR #11, fix/slicer-surface-background, 0114914 | Borrador | Procesos cancelables, optimización 100k–200k, comparación 150k/300k y evidencia de posición. Interfaz Slicer pendiente de ejecución. |
| PR #12, feat/standalone-ct-surface, 3eaf3a7 | Borrador | Script independiente por etapas, segmentación editable, proceso externo y exportación. Interfaz Slicer pendiente. |
| fix/physics-audit-20261104 | Esta revisión | Correcciones del núcleo y actualización de memoria; no fusiona ni sobrescribe las ramas anteriores. |

Las cinco PR de dependencias #4–#8 siguen abiertas. Actualizar versiones exige su propia CI; no se han mezclado con cambios físicos. La versión del paquete sigue siendo 0.5.3: el SHA identifica mejor el estado concreto.

## Flujo reconstruido

1. `rdsr.py`: lee conceptos por esquema/código y unidades UCUM, conserva valores originales, identifica eventos y comprueba totales. No usa los textos traducidos para elegir magnitudes.
2. `config.py`, `azurion.py`, `mrc200.py`: seleccionan geometría, convenciones y factores. Se registra la discrepancia de fabricante, pero se aplica la configuración escogida. El perfil Azurion conserva su control de distancias ausentes/no positivas.
3. Slicer: obtiene un **contorno corporal**, no una capa histológica de piel. Umbral, componentes y cavidades requieren inspección; una mesa conectada y brazos separados pueden invalidar el resultado. Un TC truncado genera tapas artificiales.
4. `patient.py`, `geometry.py`: carga malla, aplica unidades/ejes/origen, postura configurada y desplazamiento por evento; reconstruye fuente, detector y campo. RAS o LPS no proporcionan por sí solos el registro TC–intervención.
5. `project` y `raytrace`: evalúan centroides/puntos, divergencia, primera intersección y mesa. RTree es la ruta general; VTK es opcional para mallas densas.
6. `engine.py`: calcula por evento, aplica factores y suma sobre los mismos puntos anatómicos. Un barrido requiere aceptar explícitamente reparto uniforme; el RDSR no aporta necesariamente la secuencia temporal.
7. `metrics.py`, `reporting.py`, `interactive.py`: PSD puntual, medias locales aproximadas, áreas y exportación. La PSD es el máximo de la suma de mapas, no la suma de máximos independientes.

## Errores reproducidos y corregidos

| Gravedad | Problema anterior | Efecto y corrección |
|---|---|---|
| Alta | `field_reference` del JSON podía reinterpretar dimensiones DICOM en el isocentro. | Los conceptos 113788/113789 están definidos en el detector. Ahora el evento leído conserva `field_reference=detector`; proyección, esquema visual y comprobación KAP/área usan la misma referencia. La caché distingue referencias. Eventos manuales pueden declarar isocentro o heredar el JSON. |
| Alta en API | `calculate` aceptaba UID repetidos, aunque `read_rdsr` ya los rechazaba. | Una lista repetida podía duplicar dosis sin aviso. Se rechaza antes de cargar la malla. No se eliminan eventos silenciosamente. |
| Alta de compatibilidad | Aluminio solo reconocido como SRT C-120F9. | RDSR con SCT 12503006 o SRT C-12000 no se ejecutaban. Se admite el código actual y sus equivalentes; el proveedor espectral y el inspector Azurion comparten el mismo conjunto, evitando omitir Al añadido tras leerlo. También se admite Cu SRT C-12700. |
| Alta de presentación | Media local de 1 cm² igual a cero pese a dosis positiva. | En rejilla de 2,5 mm se descartaban los discos interiores por área y solo quedaban discos no irradiados. Ahora se devuelve métrica no disponible por muestreo; no se confunde con cero. |
| Media | Área negativa sin lados llegaba a `sqrt`. | Error genérico en lugar de diagnóstico físico. Ahora se produce InputError con UID y área no positiva. |
| Media, memoria | La suma apilaba todos los mapas en una segunda matriz. | Se conserva el orden de acumulación y los mapas por evento, usando un único vector adicional. |

Ejemplo analítico del primer error: un ancho de 200 mm en detector, SID=1000 mm, SOD=700 mm y profundidad axial de 500 mm corresponde a 100 mm de ancho. Interpretarlo en isocentro daría 142,86 mm: +42,86% en lado y +104,08% en área. No se trata de una corrección de calibración. Las configuraciones incluidas usan detector y conservan sus resultados.

No se convierte automáticamente el campo real en un campo medido distinto. Cuando solo se informa área, la hipótesis histórica de cuadrado equivalente sigue presente y ahora queda visible en cada informe de evento como `field_shape_assumption`. El área no permite recuperar aspecto, desplazamiento ni forma del campo: debe estudiarse su sensibilidad o aportar colimación adicional.

## Fundamento físico y limitaciones que permanecen

La cadena implementada es, por evento y punto visible:

`D_i(x) = Ka,r_i × (d_IRP_i / |x − foco_i|)^2 × CF_i × perfil_i(x) × T_mesa_i(x) × B_i(x) × f_i(x)`.

Los factores no utilizados en una etapa se omiten; la mesa solo se aplica si el segmento fuente–piel la cruza. La distancia de inverso del cuadrado es euclídea; el área del campo perpendicular al haz usa profundidad axial. No se añade un coseno extra a la dosis superficial: la oblicuidad de retrodispersión es una limitación física independiente, no una corrección geométrica arbitraria.

- **Origen de datos:** Ka,r, KAP, distancias y ángulos pueden proceder del RDSR; los valores faltantes normalizados desde JSON son hipótesis configuradas y llevan etiqueta. Especificaciones nominales del fabricante y espectros SpekPy no son medidas de instalación. El perfil Azurion exige distancias válidas por evento.
- **Registro:** la orientación del TC y la postura nominal RDSR no fijan el isocentro dentro del TC, la posición de brazos, los movimientos ni una deformación. La PR #11 exporta evidencia y declara que no establece registro. No se puede corregir esta carencia inventando una transformación.
- **Mesa/colchón:** transmisión escalar no modela endurecimiento espectral ni toda dispersión de mesa; verificar qué incluye la calibración para no duplicar factores. Se conservan controles contra factores incompatibles y doble BSF/conversión.
- **Benmakhlouf:** formalismo de kerma en aire incidente, retrodispersión definida en aire y conversión superficial a agua. No confundir con BSF definido como cociente de kerma en agua. La tabla se traslada localmente desde fantoma de agua; no representa anatomía, bordes, oblicuidad ni profundidad de piel individual. La memoria debe hablar de agua cuando el proveedor está en agua.
- **SpekPy:** se integra el producto monoenergético con pesos de kerma; el factor efectivo `C/B` no es el coeficiente promediado independiente de la ecuación 12. Se mantienen hashes, límites energéticos y de campo, familias de Al y tolerancia explícita de truncamiento. Una HVL compatible no determina de forma única un espectro ni valida el mapa.
- **PSD:** un centroide por cara puede perder un máximo o un borde estrecho. 150k caras no es una tolerancia dosimétrica. La media local de 1 cm² es una aproximación por grafo/vecindad, admite área muestreada distinta y no equivale automáticamente a una medida clínica. Para mallas grandes puede no calcularse; la PSD puntual sigue completa.
- **Incertidumbre:** el presupuesto parcial por cuadratura no es una incertidumbre de PSD propagada. La propagación paramétrica comparte sesgos entre eventos, pero asume independencia entre parámetros y excluye errores de anatomía y transporte. No asignar una precisión clínica a partir del número de tests.
- **Restricciones:** UID duplicado, distancias imposibles, falta de unidades y tablas fuera de dominio deben producir error. Fabricante distinto es informativo. Barrido uniforme, equivalencia cuadrada y referencias de mesa deben ser hipótesis trazables. El fallback histórico sintético, donde se activa expresamente, cambia el modelo y no valida el caso.

## Rendimiento y evidencia

La suite base pasó **163 pruebas con una omitida** en 42,55 s. Las nueve regresiones nuevas, ejecutadas contra la base, dieron **7 fallos y 2 éxitos** (dos fallos comprueban además la nueva API explícita de referencia). Tras corregir: **172 pruebas superadas y una omitida** en 39,35 s. Los tiempos de suite no son un benchmark de aceleración.

Se perfiló una carga sintética de 200 eventos, superficie analítica de 129 600 puntos y geometría repetida. El apilado consumió 0,639 s de 2,70 s con cProfile en la base. El perfil posterior no contiene ese apilado. No se atribuye toda la diferencia de ejecuciones aisladas a la corrección.

El benchmark reproducible `scripts/benchmark_accumulation.py` usa 200 mapas de 150 000 puntos, calentamiento y cinco repeticiones:

| Acumulación | Mediana | Pico adicional rastreado |
|---|---:|---:|
| Apilar y sumar | 0,1063 s | 241 201 256 bytes |
| Sumar en orden sobre vector | 0,0316 s | 1 200 200 bytes |

Diferencia máxima: **0 Gy**, igualdad exacta en esta carga y en la regresión de 31 eventos. Mejora de este paso: aproximadamente **3,36×**. Se mantienen los mapas por evento; exportar NPZ todavía los apila y puede consumir memoria. El benchmark no representa RSS total, CT real ni tiempo completo de Slicer. Resultado bruto: `benchmark_accumulation.json`.

`reproducible/check.py` conserva el mapa a rtol=1e-12 y atol=1e-15, con PSD **0,03947233203405644 Gy**, y valida hashes. `scripts/validate_benmakhlouf.py` obtiene error máximo de factor **2,22e-16** y de dosis **1,04e-17 Gy** frente al cálculo decimal independiente. Es verificación del formalismo, no validación experimental.

## Slicer y trabajo pendiente prioritario

El módulo de main ejecuta etiquetado, relleno y `CreateClosedSurfaceRepresentation` desde la llamada de interfaz. La inspección demuestra una ruta bloqueante; no demuestra que sea la causa única del cuelgue del usuario. No se ha reproducido su TC en Slicer.

Las PR #11 y #12 resuelven esa arquitectura mediante procesos externos. Se conservan como propuestas separadas: requieren ejecutar `self_test()` y un TC real en Slicer, comprobando cancelación, edición, reexportación y versiones instaladas. Una copia vieja del módulo puede seguir cargándose aunque se actualice otro checkout.

La reducción de la PR #11 evalúa integridad, desviación muestreada, área y volumen, conserva anatomía original y compara candidatos. No acredita Hausdorff continua, ausencia de todas las autointersecciones ni precisión PSD. El objetivo 100k–200k se debe aceptar solo si pasa geometría y comparación dosimétrica 150k/300k en los casos relevantes. El cierre de cavidades durante segmentación debe distinguirse de reparar agujeros para alcanzar un recuento.

Siguientes prioridades: registro medido TC–intervención; colimación asimétrica y campos solo descritos por área; transmisión/CF/HVL de instalación; validación con fantoma y película; convergencia en bordes y solapamientos; memoria de exportación y timeout del trabajador VTK. No se han inventado esos datos ni fusionado las PR antiguas para aparentar una versión validada.

## Fuentes

- DICOM PS3.16, definiciones 113788/113789/113790: https://dicom.nema.org/medical/dicom/2016d/output/chtml/part16/chapter_D.html ; detector/receptor como plano de las dimensiones.
- DICOM PS3.16, CID 10006 actual: https://dicom.nema.org/medical/dicom/current/output/chtml/part16/sect_CID_10006.html ; SCT 12503006 (Al), SCT 66925006 (Cu), equivalentes SRT.
- Benmakhlouf, Fransson y Andreo, KS-ASF-201101-IR (2011), informe aportado al proyecto: secciones 2.2, ecuaciones 1–3 y formalismo espectral. Datos originales y hashes conservados en `skindose_verifier/data/benmakhlouf_2011/`.
- Pruebas y código de las ramas indicadas; no se sustituyen medidas independientes por resultados de software.

## Comprobaciones adicionales de esta sesión

- PR #11: 21 pruebas nativas superadas (12 de malla y 9 de posición), 62,52 s, VTK 9.3.1. Se han ejecutado desde su propia rama, sin mezclar sus fuentes con main.
- PR #12: 7 pruebas nativas superadas, 3,48 s, incluido trabajador de segmentación/exportación. No son pruebas de la interfaz Qt/MRML.
- Backend VTK: 1 prueba superada con el intérprete Python de este entorno y VTK, seleccionándolo mediante `DERMADOSE_PYTHONSLICER`. Compara intersecciones axiales/oblicuas/interiores y segmentos sin intersección frente a Trimesh. No se ha ejecutado PythonSlicer ni Slicer real. En la suite ordinaria ese test se omite al no haber instalación de Slicer.
- Compilación de memoria intentada con LuaLaTeX/latexmk: no se genera PDF porque falta `spanish.ldf` (paquete de idioma español de TeX Live). Además, el repositorio no contiene varias figuras históricas referenciadas (`docs/figuras/` y `docs/evidencias/`). No se han inventado esas imágenes ni afirmado revisar un PDF inexistente. Se ha hecho portable la selección de fuente monoespaciada, manteniendo Consolas cuando esté instalada.

Para reproducir: `python -m pytest`, `python scripts/validate_benmakhlouf.py`, `python reproducible/check.py`, `python scripts/benchmark_accumulation.py`. Las pruebas de Slicer deben ejecutarse en sus ramas y después dentro de Slicer según las guías de cada PR. Restaurar las figuras históricas y el idioma español antes de compilar la memoria completa.

## Hallazgo adicional durante el estudio de convergencia

La superficie analítica sintética se evaluó a 10, 5 y 2,5 mm. Entre 5 y 2,5 mm la PSD puntual cambia +0,00795%, el área irradiada +1,545% y la integral dosis-área −0,1925%. Son magnitudes distintas: estabilidad de PSD no acredita estabilidad de la huella. Los datos se guardan en `convergencia_sintetica.json`; este estudio no es una comparación de mallas de paciente 150k/300k.

A 2,5 mm, el disco discreto interior de radio sqrt(100/pi) contiene 131,25 mm², por encima del límite de aceptación de 125 mm². Solo quedaban discos de borde no irradiados y el algoritmo devolvía 0 Gy para la media local máxima. La corrección conserva la definición existente y devuelve `None`/`unavailable_spatial_sampling` cuando todos los discos admisibles tienen cero pero el mapa contiene dosis positiva. Se aplica también a la rama de mallas. Un mapa genuinamente nulo conserva su resultado cero cuando tiene vecindades admisibles. No se sustituye esta métrica por un promedio arbitrario ni se pretende haber implementado integración exacta de 1 cm².

Se añadieron dos pruebas que reproducen esta falsa media nula en la base y comprueban la señalización de indisponibilidad en el informe corregido. Pendiente: integración parcial de celdas o refinamiento local para estimar áreas exactas, con una definición espacial consistente y verificable.

## Cierre y repetición final — 6 de noviembre de 2026

Tras reanudar la sesión se verificó que main y la rama de auditoría seguían en c59ed26, sin cambios ajenos que sobrescribir. Con las once regresiones finales, la suite completa terminó con **174 pruebas superadas y una omitida en 25,87 s**. Se repitió además `reproducible/check.py`: mapa y PSD conservados. El test omitido requiere una instalación detectable de PythonSlicer; su comprobación nativa separada con VTK se documentó arriba y no equivale a probar Slicer. Los tiempos de diferentes sesiones no se usan como medida de aceleración.

Se publican las correcciones y esta memoria en una propuesta hacia main; no se integran automáticamente las PR de Slicer ni se cambia main. El entorno final queda en `entorno.json` y la salida de pruebas en `pruebas_finales.txt`.
