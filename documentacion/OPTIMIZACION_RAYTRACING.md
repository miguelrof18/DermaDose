# Optimización del trazado de rayos

## Cambio aplicado

El detector de primera intersección usa el árbol R-tree de Trimesh. Antes, cada consulta pasaba por una lista temporal en `CountingTree.intersection` y Trimesh volvía a materializarla. Además, la búsqueda consideraba la caja de la malla a lo largo de toda la semirrecta, incluso más allá del punto de piel que se estaba evaluando.

Ahora el contador envuelve el iterador de candidatos sin crear una lista intermedia y limita la consulta al segmento entre la fuente y el punto de piel. Se conserva el lote acotado de 64 rayos para mantener acotado el uso de memoria.

## Comparación reproducible

Se procesó el mismo RDSR externo de cinco eventos, con el mismo JSON XCAT50 y factores sintéticos históricos. En ambas mediciones se usó `cProfile`:

| Medida | Antes | Después |
|---|---:|---:|
| Tiempo de `calculate` | 28,23 s | 2,65 s |
| Rayos | 5.192 | 5.192 |
| Pares rayo-triángulo candidatos | 31.675.064 | 55.357 |
| Diferencia máxima del mapa de dosis | — | 0 Gy |

El mapa calculado fue idéntico punto por punto (`numpy.array_equal`). El tiempo incluye instrumentación de perfilado; la exportación de mapas e informes no forma parte de `calculate`.

## Separación del código

El núcleo ya está dividido por responsabilidad: `rdsr.py` lee el DICOM, `geometry.py` modela el haz y proyecta la geometría, `raytrace.py` resuelve la primera intersección, `engine.py` aplica factores y acumula eventos, y `metrics.py` calcula métricas. Se mantiene esa división y se evita trasladar funciones entre módulos durante una optimización de rendimiento, para reducir el riesgo de alterar la dosimetría.

## Empaquetado

El programa genera `skindose_verifier/data/mapa_3d_interactivo.html` como archivo de salida histórico, pero el generador solo necesita `visor_3d_template.html`. El manifiesto de paquete ahora incluye únicamente esa plantilla, evitando distribuir el HTML generado de aproximadamente 34 MB junto al paquete.
