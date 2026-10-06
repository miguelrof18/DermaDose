# ROADMAP DEL CALCULADOR DE DOSIS EN PIEL

| Funcionalidad | Estado actual | Mejora propuesta | Evidencia bibliografica | Prioridad | Complejidad | Implementado | Validacion necesaria |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Lectura RDSR DICOM | Parser funcional con campos sinteticos y comprobaciones de presencia | Ampliar inventario multivendor con RDSR anonimizados reales | TG-357 [R1], DICOM PS3.17 [R15] | Alta | Media | Parcial | Revision manual de eventos por equipo y comparacion con consola/RDSR |
| Geometria fuente-isocentro-detector | Reconstruccion determinista y trazable | Validar signos y referencias para equipos reales biplano | TG-357 [R1], Johnson [R3] | Alta | Alta | Parcial | Fantoma geometrico, campos conocidos y QA por angulos |
| Malla triangular 3D | XCAT y PySkinDose disponibles; refinamiento conformal | Segmentacion de paciente real desde CT/CBCT o superficie optica | Johnson [R3], Krajinovic [R4], XCAT [R16] | Alta | Alta | Si, para fantomas | Comparar superficie real, landmarks y errores de registro |
| Refinamiento de malla | Aristas maximas 40/20/10/5 mm y ROI opcional | Refinamiento adaptativo guiado por borde de campo y gradiente | Johnson [R3], Trimesh [R17] | Alta | Media | Si, primera version | Convergencia en geometrias analiticas y fantomas medidos |
| Ray tracing | Primer impacto con Rtree y contadores | Evaluar Embree/Open3D/VTK si el tamano de malla crece mucho | Rana [R6], Takata [R7] | Media | Media/Alta | Si, CPU actual | Benchmark con mallas >1M caras y mapas equivalentes |
| Factores dosimetricos | Interpolacion JSON y ejes de angulo/campo | Sustituir constantes sinteticas por curvas medidas de cada equipo | Ganan Mora [R2], Malchair [R5], IAEA TRS-457 [R19] | Alta | Media | Parcial | Calibracion experimental, incertidumbre y rango de validez |
| Angulo de incidencia | Variable calculada y disponible para factores | Medir o modelar dependencia angular real por calidad de haz | Ganan Mora [R2], Omar [R11] | Alta | Media | Parcial | Medidas en fantoma inclinado y comparacion con literatura |
| Backscatter y conversion a piel | Factores configurables, no clinicos | Tablas dependientes de energia, campo, material y geometria | IAEA TRS-457 [R19], ICRU 74 [R20] | Alta | Media | Parcial | Medidas con camara/film y presupuesto de incertidumbre |
| Scatter interno y heterogeneidades | No modelado | Mantener fuera del nucleo o comparar con Monte Carlo externo | EGSnrc [R22], TOPAS/GATE [R23], PENELOPE [R24] | Media | Alta | No | Comparacion MC en casos canonicos |
| Movimiento/deformacion paciente | No modelado | Registrar mesa, respiracion o superficie durante el procedimiento | Krajinovic [R4], future work revisado | Media | Alta | No | Datos temporales reales y prueba de acumulacion deformable |
| Metricas PSD/Dmean/D2/D5/D10 | Implementadas con ponderacion por area | Definir conjunto de metricas clinico-investigador por protocolo | Malchair [R5], TG-357 [R1] | Alta | Baja | Si | Concordancia con software de referencia y definiciones aprobadas |
| Histogramas dosis-area | CSV y PNG generados | Exportar comparadores normalizados y curvas por evento | Malchair [R5] | Media | Baja | Si | Test de monotonia y comparacion con mapas medidos |
| Hotspots | Clustering por conectividad de caras >=90 % | Umbrales configurables y seguimiento temporal por evento | Johnson [R3], Krajinovic [R4] | Alta | Baja/Media | Si | Fantomas con dos zonas calientes conocidas |
| Visor 3D | HTML offline con cortes, fuente, haz, gantry e inspeccion | Verificacion visual WebGL en Windows y modo ligero para mallas grandes | Rana [R6], tesis PyMCGPU-IR [R9] | Alta | Media | Si | Capturas desktop/movil, prueba de no solapamiento y carga real |
| Incertidumbre | LHS/MC parametrico con semilla | Presupuesto completo con correlaciones y componentes experimentales | JCGM 101 [R14] | Alta | Alta | Parcial | Datos de calibracion, repetibilidad y validacion estadistica |
| Validacion automatica | 52 tests de geometria, fisica sintetica y regresion | Anadir suite con RDSR reales anonimizados y golden files | TG-357 [R1], metodologia RI [R2] | Alta | Media | Si | CI reproducible y revision de tolerancias |
| Comparacion con TPS/MC | Documentada como protocolo, no ejecutada | Ejecutar EGSnrc/TOPAS/PENELOPE o TPS con fantomas equivalentes | EGSnrc [R22], TOPAS/GATE [R23], PENELOPE [R24] | Media | Alta | No | Comparacion independiente de mapas y PSD |
| Validacion experimental | Protocolo escrito | Medir campo, huella, PSD y area con film/TLD/OSLD/camara | Ganan Mora [R2], IAEA TRS-457 [R19] | Critica | Alta | No | Campana dosimetrica por equipo y tolerancias |
| Reproducibilidad | Scripts, hashes, versiones y lock informativo | Empaquetado con dataset de prueba versionado y CI | JCGM [R14], buenas practicas software | Alta | Baja/Media | Si | Reejecucion en Windows del usuario y entorno limpio |
| Modelos surrogate/ML/GNN | Solo revisado como idea de investigacion | Usar solo cuando haya dataset medido/MC suficiente | Han et al. 2025 [R12], lineas futuras [R9, R10] | Baja/Media | Muy alta | No | Dataset externo, validacion fuera de muestra y explicabilidad |
| GPU ray tracing | No implementado | Evaluar si el calculo supera limites CPU en casos reales | Rana [R6], Takata [R7] | Media | Alta | No | Benchmark CPU/GPU con igualdad numerica aceptable |
| SDF/differentiable rendering | Idea futura | Aplicar a registro anatomico o optimizacion, no al nucleo actual | Revision de metodos mesh/SDF [R17] | Baja | Muy alta | No | Pruebas de exactitud geometrica y estabilidad |

## Criterio de avance

La siguiente fase deberia priorizar validacion experimental y datos reales antes de anadir complejidad fisica. El nucleo actual ya permite investigar geometria, mallado, sensibilidad y visualizacion; lo que falta para uso clinico no es solo codigo, sino calibracion, tolerancias, QA y comparacion independiente.
