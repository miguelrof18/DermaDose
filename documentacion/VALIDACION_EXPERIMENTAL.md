# Protocolo propuesto de validación externa

**No se ha ejecutado una validación experimental en esta fase.** Se han verificado implementación, invariantes geométricos, regresiones y sensibilidad numérica. Las medidas de las publicaciones no se han reutilizado como si procedieran de nuestro programa.

## Magnitud y registro antes de comparar

Especificar kerma en aire libre, kerma incidente o dosis absorbida, medio, profundidad de evaluación y huella del dosímetro. Una dosis de entrada aproximada sobre un mesh sin espesor no es equivalente automáticamente a dosis a 0,07 mm, a una capa de 2 mm ni a la dosis de un vóxel. Documentar calibración de Ka,r, kVp/HVL/filtración, campo en su plano de referencia, mesa/colchoneta y sistema de coordenadas. Registrar el fantoma en milímetros al mismo marco que `patient_points_mm`.

## Secuencia experimental

| Bloque | Ensayo | Referencia y decisión que permite |
| --- | --- | --- |
| Metrología de partida | Cámara trazable en aire: Ka,r, linealidad con exposición y distancia. | Distingue calibración de errores geométricos; rama `tube_output` independiente cuando su rendimiento es medido. |
| Mesa | Mismo haz con/sin mesa y colchoneta; diferentes kVp/Cu y ángulos. | Separar transmisión efectiva y su dependencia; evitar doble conteo de BSF/scatter. |
| Plano de agua/equivalente | Distancias, campos, ángulos y traslaciones, primero por separado y luego combinados. | Esquema de Gañán Mora et al. [R2]. La incertidumbre del 13% (k=2) de su trabajo NO se asigna a nuestras medidas. |
| Perfiles y bordes | Película radiocrómica calibrada para la calidad del haz, con registro y curva incertidumbre. | Perfil/heel, penumbra, superficie de isodosis, contribución fuera del campo. |
| Curvatura | Cilindro y fantoma antropomórfico; vistas PA/AP/laterales/oblicuas; brazos. | Primera entrada, normales, autosombra y correspondencia anatómica. |
| Acumulación | Campos coincidentes, parcialmente solapados y separados; cambios de mesa. | La suma de máximos individuales no es la PSD acumulada salvo coincidencia espacial. |
| Repetibilidad | Montajes repetidos, perturbaciones conocidas y diferentes operadores. | Estimar distribuciones/correlaciones para el módulo de incertidumbre. |
| Compatibilidad | RDSR reales anonimizados de cada fabricante/modelo/versión. | Conformidad de códigos/unidades y signos, origen de mesa, tamaño de campo y filtros. |

Las cámaras proporcionan una referencia puntual con volumen sensible finito; la medida puede requerir correcciones de perturbación y conversión. Película permite distribución bidimensional, pero necesita calibración energética, lectura reproducible y registro. TLD/OSLD permiten muestreo distribuido; considerar respuesta angular/energética, área sensible y localización. No convertir una coincidencia numérica en aceptación sin un presupuesto metrológico.

## Herramientas ya disponibles

`validate-measurements` recibe `test_id,rdsr_path,config_path,measured_max_gy,standard_uncertainty_gy` y devuelve diferencias relativas y normalizadas por la incertidumbre de la medida. Esta última no incluye la del modelo. `validate-map` recibe puntos registrados `x_mm,y_mm,z_mm,measured_gy,standard_uncertainty_gy,quantity`, calcula RMSE, separación de picos y Dice en los puntos muestreados; exige distancia máxima y asociaciones no duplicadas. No calcula gamma ni Dice de toda la superficie. Los CSV vacíos de `examples/` son plantillas, sin observaciones inventadas.

## Monte Carlo, TPS y otros verificadores

Geant4/TOPAS, EGSnrc o PENELOPE pueden actuar como referencias externas [R22–R24]. Antes de comparar, modelar el espectro/normalización, mesa y fantoma; fijar materiales, espesor de piel, scorer, cortes de producción y estadística. Comparar primero dosis/kerma en geometrías simples; después registrar la superficie común y estudiar mapas, áreas y picos. La PSD obtenida como máximo de muchas celdas MC puede sesgarse por ruido [R12]; no comparar contra un máximo ruidoso sin estudiar su estabilidad.

GATE ofrece infraestructura Geant4. PRIMO y TPS de radioterapia pueden servir para preguntas MV específicas, pero no se asume que modelen un arco de intervencionismo de kV. Los DICOM RTPLAN/RTDOSE no se leen en este programa. `compare-radimetrics` conserva la comparación comercial como secundaria: dos programas que comparten datos o hipótesis pueden coincidir y estar ambos sesgados.

Las tolerancias de aceptación deberán fijarse para cada magnitud, rango y uso previsto, incorporando la incertidumbre de medida, registro, discretización y modelo. No se heredan ±25%, 13% ni otros umbrales de artículos como una tolerancia universal del calculador.
