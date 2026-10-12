# CHANGELOG_RESEARCH.md

Trazabilidad de la fase de I+D v0.5.0. Las referencias `[R...]` se describen en `RESEARCH_REVIEW.md`.

## 2026-09-28 - v0.5.0

### Auditoria tecnica y bibliografia

- **Descripcion:** se reviso el modelo fisico, geometrico, numerico, de rendimiento y de visualizacion.
- **Motivacion:** el prompt solicitaba una fase de I+D, no una lista teorica de mejoras.
- **Referencias:** AAPM TG-357 [R1], Johnson et al. [R3], Krajinovic et al. [R4], Malchair et al. [R5], Rana et al. [R6], Takata et al. [R7], Omar et al. [R8], JCGM 101 [R14], DICOM PS3.17 [R15].
- **Archivos:** `AUDITORIA_TECNICA.md`, `RESEARCH_REVIEW.md`, `VALIDACION_EXPERIMENTAL.md`.
- **Validacion:** se documento cada hallazgo con impacto, solucion, prioridad y prueba asociada.

### Refinamiento de malla triangular

- **Descripcion:** nuevo refinamiento conformal por aristas compartidas con conservacion de area/volumen y opcion de ROI.
- **Motivacion:** el parametro de resolucion anterior describia la malla, pero no la refinaba realmente.
- **Referencias:** literatura de skin-dose mapping sobre dependencia geometrica [R3, R4, R6, R7] y documentacion Trimesh [R17].
- **Archivos:** `skindose_verifier/mesh_processing.py`, `patient.py`, `convergence.py`, `cli.py`.
- **Efecto esperado:** reducir sensibilidad al tamano irregular de triangulos en bordes de campo y hotspots.
- **Validacion:** pruebas de conservacion, winding, ROI, limite de caras y convergencia XCAT 40/20/10/5 mm.

### Ray tracing instrumentado

- **Descripcion:** se aislo el primer impacto rayo-triangulo y se anadieron contadores de rayos y pares candidato.
- **Motivacion:** medir rendimiento real y evitar dobles impactos por superficies ocultas.
- **Referencias:** metodos de interseccion superficie-haz en skin-dose mapping [R3, R6, R7].
- **Archivos:** `skindose_verifier/raytrace.py`, `geometry.py`, `engine.py`.
- **Efecto esperado:** resultados mas trazables y benchmark interpretable.
- **Validacion:** pruebas de primera interseccion, autosombra, cilindro inclinado, simetrias y contadores.

### Correcciones geometricas del haz

- **Descripcion:** calculo de area de campo a profundidad axial, angulo de incidencia local, coordenadas normalizadas de haz y rechazo de fuentes dentro de la malla.
- **Motivacion:** evitar sesgos en campos oblicuos y geometrias degeneradas.
- **Referencias:** DICOM PS3.17 para marcos de coordenadas [R15], TG-357 para trazabilidad geometrica [R1].
- **Archivos:** `skindose_verifier/geometry.py`, `config.py`, `engine.py`.
- **Efecto esperado:** campo proyectado y factores dependientes de angulo/campo mas coherentes.
- **Validacion:** tests de area off-axis, signos, reflexion por angulo opuesto, cilindro inclinado y rotacion de malla.

### Factores dosimetricos configurables

- **Descripcion:** interpolacion cacheada de factores y nuevos ejes `skin_angle_deg`, `field_fraction_x` y `field_fraction_y`.
- **Motivacion:** permitir calibraciones experimentales sin codificar constantes ficticias.
- **Referencias:** Ganan Mora et al. [R2], Malchair et al. [R5], IAEA TRS-457 [R19].
- **Archivos:** `skindose_verifier/config.py`, `engine.py`, `data/equipo_sintetico.json`.
- **Efecto esperado:** base preparada para curvas medidas de BSF, conversion y perfil de haz.
- **Validacion:** pruebas analiticas con formulas cerradas sinteticas.

### Metricas avanzadas de dosis

- **Descripcion:** histogramas dosis-area, Dmean, Dmax, D2%, D5%, D10%, areas por umbral, integral dosis-area y hotspots conectados.
- **Motivacion:** la PSD puntual no describe extension espacial ni regiones de alta dosis.
- **Referencias:** revision de software de dosis cutanea [R5], aplicaciones independientes [R4], JCGM para trazabilidad de incertidumbre [R14].
- **Archivos:** `skindose_verifier/metrics.py`, `engine.py`, `reporting.py`.
- **Efecto esperado:** salida mas util para investigacion, comparacion y QA.
- **Validacion:** tests de cuantiles ponderados por area, histogramas monotono y hotspots separados.

### Incertidumbre parametrica

- **Descripcion:** modulo `uncertainty` con Monte Carlo ordinario o Latin Hypercube Sampling, semillas y salida JSON/NPZ.
- **Motivacion:** separar mapa determinista de sensibilidad a posicion, angulos y factores dosimetricos.
- **Referencias:** JCGM 101 [R14], SciPy LHS [R18], validacion RI [R2].
- **Archivos:** `skindose_verifier/uncertainty.py`, `cli.py`, `examples/incertidumbre_sintetica.json`.
- **Efecto esperado:** primer mapa probabilistico y distribucion de PSD bajo hipotesis explicitas.
- **Validacion:** tests de semilla, escala lognormal, incertidumbre cero y desplazamiento de paciente; ejecuciones de 64 y 256 muestras.

### Benchmark y memoria

- **Descripcion:** benchmark con/sin cache geometrica, escalabilidad por caras, memoria RSS en procesos limpios y contadores de ray tracing.
- **Motivacion:** cuantificar cuellos de botella antes de considerar GPU/Embree/Numba.
- **Referencias:** Rana et al. GPU [R6], Takata et al. [R7], tesis PyMCGPU-IR [R9, R10].
- **Archivos:** `skindose_verifier/benchmark.py`, `scripts/measure_memory.py`, `scripts/reproduce_research.py`.
- **Efecto esperado:** cache geometrica reutiliza proyecciones identicas sin cambiar el mapa.
- **Validacion:** benchmark XCAT 10 mm: 8,94x en caso de 12 campos identicos, diferencia maxima 0,0 Gy.

### Visor 3D interactivo

- **Descripcion:** HTML Plotly offline con cortes, eventos, fuentes, isocentro, detector, rayos, campo, arco gantry, hotspots, DAH e inspeccion de caras.
- **Motivacion:** el usuario pidio poder ver el mapa 3D moviendo planos y perspectivas.
- **Referencias:** Rana et al. [R6], Johnson et al. [R3], tesis PyMCGPU-IR [R9].
- **Archivos:** `skindose_verifier/interactive.py`, `skindose_verifier/data/visor_3d_template.html`, `scripts/check_viewer.cjs`.
- **Efecto esperado:** exploracion espacial mucho mas clara sin alterar los datos numericos.
- **Validacion:** prueba DOM/Plotly stub superada sobre HTML de 177 926 caras. Falta inspeccion WebGL real en navegador.

### Reproducibilidad cientifica

- **Descripcion:** runtime audit con dependencias y hashes, scripts de reproduccion, lock informativo, plantillas CSV y documentacion cientifica.
- **Motivacion:** que los resultados se puedan repetir y auditar.
- **Referencias:** buenas practicas de software cientifico y documentos de validacion [R1, R2, R14].
- **Archivos:** `provenance.py`, `requirements-research-lock.txt`, `scripts/reproduce_research.py`, `README.md`, `INICIO_AQUI.txt`.
- **Efecto esperado:** ejecuciones con version, semilla, fuente de datos y limitaciones visibles.
- **Validacion:** logs reproducibles en `research_results/` y 52 pruebas automatizadas.

## v0.5.1 - cambios derivados de notas de revisión

### Localización automática de RDSR

- **Descripcion:** nuevo comando `find-rdsr` para buscar informes `X-Ray Radiation Dose SR` dentro de una carpeta DICOM y diferenciar imágenes XA/RF de informes de dosis.
- **Motivacion:** facilitar el uso con carpetas clínicas reales, donde el procedimiento contiene múltiples series y no siempre es evidente qué fichero debe entrar al calculador.
- **Referencias:** DICOM PS3.16/PS3.17 para RDSR e Irradiation Event X-Ray Data [R15].
- **Archivos:** `skindose_verifier/rdsr.py`, `skindose_verifier/cli.py`, `tests/test_pipeline.py`, `README.md`.
- **Efecto esperado:** flujo más robusto desde exportaciones DICOM reales sin leer píxeles ni exponer PatientName/PatientID.
- **Validacion:** prueba automatizada con un RDSR sintético y un DICOM XA no candidato; ejecución manual sobre `results/demo_xcat`.

### Contribución de eventos al punto de PSD

- **Descripcion:** `resultados.json` incorpora `global_psd_event_contributions`, ranking de eventos según la dosis aportada en el punto de máxima dosis acumulada.
- **Motivacion:** convertir la información evento a evento en interpretación dosimétrica, no solo en suma final.
- **Referencias:** mapas de dosis por evento y análisis de acumulación en software de dosis cutánea [R3, R4, R5].
- **Archivos:** `skindose_verifier/engine.py`, `tests/test_pipeline.py`, `README.md`.
- **Efecto esperado:** identificar si el máximo procede de un evento dominante o de superposición de campos.
- **Validacion:** prueba automatizada que verifica ranking y suma exacta de contribuciones en el punto de PSD.

## v0.5.2 - trazabilidad dosimétrica de factores

### Metadatos completos de factores

- **Descripcion:** `resultados.json` incorpora `factor_metadata`, con procedencia, identificador, representación, valor o ejes de interpolación, e incertidumbre relativa estándar de cada factor.
- **Motivacion:** responder a las notas de revisión que pedían distinguir factores sintéticos, medidos, tabulados o calculados, y documentar dominio e incertidumbre.
- **Referencias:** Benmakhlouf et al. para factores tabulados de retrodispersión/conversión, TG-357 y protocolos de validación para trazabilidad de factores [R1, R2, R25].
- **Archivos:** `skindose_verifier/engine.py`, `tests/test_pipeline.py`, `README.md`.
- **Efecto esperado:** preparar el calculador para incorporar tablas reales de BSF, transmisión de mesa y conversión aire--piel sin confundirlas con valores sintéticos.
- **Validacion:** prueba automatizada sobre presencia, procedencia, identificador y tipo de representación de los factores.

## v0.5.3 - informe de compatibilidad y trazabilidad de notas

### Informe de compatibilidad RDSR

- **Descripcion:** nuevo comando `compatibility-report` para revisar candidatos RDSR antes de ejecutar el cálculo completo.
- **Motivacion:** responder a la nota sobre aceptar RDSR reales anonimizados y generar un informe de compatibilidad.
- **Archivos:** `skindose_verifier/cli.py`, `tests/test_pipeline.py`, `README.md`.
- **Efecto esperado:** identificar si el parser lee el RDSR, si faltan campos críticos, si falta calidad de haz para la etapa solicitada y si el fabricante/modelo coincide con la configuración.
- **Validacion:** prueba automatizada con RDSR sintético y configuración coincidente.

### Trazabilidad de comentarios `% Nota:`

- **Descripcion:** nuevo documento `NOTAS_IMPLEMENTACION.md` que agrupa las anotaciones `% Nota:` por tema y declara estado de implementación.
- **Motivacion:** dejar claro qué se implementó, qué quedó preparado y qué requiere calibración o datos reales.
- **Archivos:** `NOTAS_IMPLEMENTACION.md`.
- **Efecto esperado:** evitar que comentarios ocultos en LaTeX queden sin respuesta técnica.

## Pendiente tras v0.5.3

- Validacion experimental con fantoma, peliculas/TLD/OSLD/camara y RDSR reales anonimizados.
- Comparacion externa con TPS, Monte Carlo o calculadores comerciales.
- Curvas dosimetricas medidas por equipo, modo, energia, filtracion, campo, mesa y colchoneta.
- Evaluacion WebGL real en Windows/VS Code del usuario.
- Segmentacion de pacientes reales y control de auto-intersecciones.
- Presupuesto de incertidumbre completo con correlaciones y componentes de calibracion.
