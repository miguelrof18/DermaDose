# Acta de verificacion computacional v0.5.0

Fecha de actualizacion: 2026-09-28. Proyecto: `skindose-verifier` 0.5.0.

Este documento resume la verificacion computacional de la fase de I+D. No constituye validacion fisica ni clinica. Los coeficientes del equipo, el RDSR y los escenarios de incertidumbre son sinteticos salvo que se indique lo contrario.

## Comprobaciones completadas

- `python3 -m pytest -q`: **52 pruebas superadas**.
- Se cubren: lectura RDSR, geometria del haz, campo proyectado, ley inversa del cuadrado, factores por modo/calidad/campo/angulo, malla triangular, cilindro, esfera, autosombra, signos de angulos, acumulacion, mapas, histogramas dosis-area, hotspots, incertidumbre LHS, remallado, limites de entrada y errores por NaN/Inf.
- El estudio reproducible `scripts/reproduce_research.py --skip-tests` finalizo correctamente y genero los informes en `research_results/`.
- La logica del visor HTML se comprobo con `scripts/check_viewer.cjs`: controles, cortes X/Y/Z, datos de eventos, inspeccion de puntos, histograma, hotspots y foco regional. Esta prueba no sustituye una inspeccion WebGL real en navegador.
- Se construyo y ejecuto el ejemplo XCAT refinado con malla triangular de piel procedente del repositorio solicitado, usando arista maxima de 10 mm.

## Resultado principal XCAT refinado

| Indicador | Resultado |
| --- | ---: |
| Fantoma | XCAT WholeBody Body STL |
| Modelo | Malla triangular 3D cerrada |
| Arista maxima nominal | 10 mm |
| Caras de superficie | 177 926 |
| Eventos sinteticos | 2 |
| PSD puntual | **0,034776306165 Gy** |
| Media local tipo 1 cm2 | 0,034326648816 Gy |
| Area real de esa vecindad | 89,851828823 mm2 |
| Superficie irradiada | 139,855228184 cm2 |
| Hotspots >= 90 % del maximo | 1 region |
| Area hotspot principal | 104,789097179 cm2 |
| Centro del maximo | (-51,257; -30,033; -96,623) mm |
| Rayos trazados, 2 eventos | 1 959 |
| Pares candidato rayo-triangulo | 216 970 |

La media local no se fuerza a exactamente 1 cm2: se informa el area real agregada porque depende del grafo de caras y de la discretizacion local. El resultado demuestra coherencia numerica sobre una malla realista, no exactitud dosimetrica frente a un paciente o equipo real.

## Convergencia de malla

Estudio realizado sobre la misma anatomia, refinando aristas compartidas sin cambiar la forma macroscopica del STL.

| Arista maxima | Caras | PSD puntual (Gy) | Area irradiada (cm2) | Integral dosis-area (Gy cm2) | Tiempo (s) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 40 mm | 16 390 | 0,034720081694 | 142,447935585 | 4,151890628 | 0,188 |
| 20 mm | 46 886 | 0,034761885529 | 141,417302758 | 4,220602277 | 0,525 |
| 10 mm | 177 926 | 0,034776306165 | 139,855228184 | 4,185788510 | 2,234 |
| 5 mm | 709 122 | 0,034782052219 | 140,408175140 | 4,200290426 | 10,547 |

Cambio 10 mm -> 5 mm:

- PSD puntual: +0,0165 %.
- Area irradiada: +0,395 %.
- Integral dosis-area: +0,346 %.
- Desplazamiento del punto maximo: 2,333 mm.

Interpretacion: la PSD puntual es estable en este caso sintetico. Las magnitudes de area cambian mas porque dependen de los bordes del campo y del mallado; no se asume convergencia monotona.

## Rendimiento y memoria

Benchmark sobre 12 campos identicos, superficie ya preparada y resumen minimo:

| Caso | Mediana sin cache | Mediana con cache | Aceleracion | Diferencia maxima de mapa |
| --- | ---: | ---: | ---: | ---: |
| XCAT 10 mm, 177 926 caras | 1,311 s | 0,147 s | 8,94x | 0,0 Gy |

Medida RSS en procesos limpios, excluyendo exportacion HTML/PNG:

| Arista maxima | Caras | Pico RSS | Tiempo calculo |
| ---: | ---: | ---: | ---: |
| Original | 12 148 | 95,4 MB | 0,083 s |
| 20 mm | 46 886 | 124,4 MB | 0,553 s |
| 10 mm | 177 926 | 227,2 MB | 2,430 s |
| 5 mm | 709 122 | 642,8 MB | 11,622 s |

Estas cifras son dependientes del equipo de ejecucion y del caso sintetico. Sirven para orientar escalabilidad, no para prometer rendimiento universal.

## Incertidumbre parametrica

Primer modulo de incertidumbre por Latin Hypercube Sampling, con 64 muestras y semilla `20260928`. Parametros hipoteticos: desplazamiento paciente 3 mm por eje, sesgo angular 0,5 grados, escala Ka lognormal 5 % y BSF lognormal 3 %.

| Indicador PSD | Resultado |
| --- | ---: |
| Media | 0,034789350163 Gy |
| Desviacion tipica | 0,002088140176 Gy |
| p2,5 | 0,031431119134 Gy |
| Mediana | 0,034718160455 Gy |
| p97,5 | 0,039526007931 Gy |

Una ejecucion adicional de 256 muestras con otra semilla dio media 0,034782578479 Gy y desviacion tipica 0,002107513269 Gy. Los cuantiles siguen siendo sensibles al muestreo y a las hipotesis; no son un presupuesto de incertidumbre clinico.

## Visor 3D

El visor `results/research_xcat/mapa_3d_interactivo.html` incluye:

- mapa de dosis acumulada y por evento;
- cortes X/Y/Z ajustables;
- escala fija o automatica;
- fuente, isocentro, detector, campo, rayos principales y arco esquematico del gantry;
- hotspots y coordenadas;
- histograma dosis-area acumulado;
- area que supera un umbral elegido;
- inspeccion de cara: coordenadas, normal, dosis total, dosis por evento, distancia e incidencia cuando hay geometria disponible;
- opcion de color por vertices solo visual, sin modificar los valores numericos.

No pude ejecutar una prueba grafica real con Chromium/Firefox en este entorno porque no habia navegador disponible y la descarga de Playwright fallo. La prueba automatizada cubre la logica DOM/Plotly, pero la inspeccion WebGL final debe hacerse abriendo el HTML en el ordenador del usuario.

## Limites actuales

- No hay validacion con peliculas, TLD/OSLD, camara, TPS o Monte Carlo.
- No hay calibracion real de equipo ni curvas reales de conversion, BSF, mesa o colchoneta.
- El modelo calcula dosis de entrada/superficie con aproximaciones deterministas; no transporta fotones/electrones en tejido.
- La malla es un fantoma XCAT, no una segmentacion personalizada del paciente.
- No se corrigen heterogeneidades internas, organos, movimiento respiratorio ni deformacion de paciente.
- La deteccion de auto-intersecciones del STL no se resuelve automaticamente.

## Conclusion

La version 0.5.0 cumple la fase computacional pedida: malla triangular 3D, refinamiento real, estudio de convergencia, ray tracing instrumentado, cache geometrica, metricas avanzadas de dosis en piel, visor interactivo ampliado, incertidumbre parametrica y 52 pruebas automatizadas. Sigue siendo software de investigacion y verificacion independiente; para uso clinico requiere validacion dosimetrica especifica, comparacion con metodos de referencia, tolerancias y QA.
