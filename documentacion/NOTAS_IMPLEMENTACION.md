# Trazabilidad de comentarios `% Nota:`

Este documento resume las anotaciones del usuario que empiezan exactamente por
`% Nota:` en la memoria modular aportada y su estado en el calculador.

## Resumen inicial

| Nota | Estado | Implementación / decisión |
| --- | --- | --- |
| Mejorar visor con escala configurable, umbral `A(D>=D0)` y curvas dosis-area | Implementado parcialmente | El visor y las salidas ya incluyen histograma dosis-area, CSV de curva acumulada, mapa interactivo y exploración por evento. La escala fija/relativa queda documentada; no se recalcula dosis al cambiar umbral. |
| Mejorar métricas: areas por umbrales, contribución por evento, ranking de eventos | Implementado | `surface_metrics`, `histograma_dosis_area.csv`, `event_dose_gy` y `global_psd_event_contributions`. |
| Mejorar trazabilidad de factores | Implementado | `factor_metadata` informa procedencia, identificador, tipo, valor/ejes e incertidumbre relativa. |
| Buscar automáticamente el RDSR en carpeta DICOM | Implementado | Comando `find-rdsr`. |
| Aceptar RDSR reales anonimizados y generar informe de compatibilidad | Implementado como pre-chequeo | Comando `compatibility-report`, sin ejecutar cálculo de dosis completo ni validar físicamente. |

## Métricas dosis-area y visor

| Nota | Estado | Implementación / decisión |
| --- | --- | --- |
| Area total irradiada y areas que superan umbrales | Implementado | `surface_dose_area_approx`, `surface_metrics` e histograma dosis-area. |
| Sumar areas de caras por encima de umbral | Implementado | El cálculo usa `area_mm2` por elemento superficial. |
| Métricas por evento y acumuladas | Implementado | Mapas por evento en `event_dose_gy`, resumen por evento y ranking en el punto PSD. |
| Mapas, histogramas y curvas acumuladas | Implementado | `mapa_2d.png`, `mapa_3d.png`, `mapa_3d_interactivo.html`, `histograma_dosis_area.csv/png`. |
| Ejemplo numérico 1,2 Gy / 35 cm2 / 8 cm2 | No incorporado como resultado | Es un ejemplo ilustrativo sin caso real asociado; no se introduce como resultado del programa. |
| Barra deslizante e isolíneas en visor | Parcial | El visor interactivo permite exploración espacial/eventos; la curva dosis-area se exporta. Isolíneas y barra de umbral quedan como mejora gráfica futura. |
| Exportación de capturas o tablas de `A(D>=D0)` | Implementado para tablas | CSV y PNG de histograma dosis-area. Capturas manuales del visor quedan fuera del motor. |

## Monte Carlo de parámetros

| Nota | Estado | Implementación / decisión |
| --- | --- | --- |
| Monte Carlo de parámetros, no transporte | Implementado | Módulo `uncertainty`, LHS/Monte Carlo ordinario, semillas y mapas probabilísticos. |
| Variar dosis máxima, area irradiada y posición del máximo | Implementado parcialmente | Se calculan PSD por muestra, cuantiles, mapas medios/desviación/probabilidad; sensibilidad de posición disponible con `sensitivity-position`. |
| No sustituye validación experimental | Documentado | Marcado como análisis de incertidumbre parcial, sin validez clínica. |

## Calibración mesa-equipo

| Nota | Estado | Implementación / decisión |
| --- | --- | --- |
| Calibración geométrica mesa-equipo | Preparado, no ejecutado | La geometría de mesa se define explícitamente en configuración; requiere medidas locales. |
| Transmisión por ángulo, calidad y zona atravesada | Preparado | Factores tabulados permiten ejes como `kvp`, `filtration_mm_cu`, `table_angle_deg`, `field_fraction_x/y`. |
| Registrar fecha, método, incertidumbre y rango | Implementado a nivel de metadatos | `factor_metadata` conserva procedencia, identificador, incertidumbre y ejes/rangos. Fecha/método deben codificarse en `identifier`/`provenance` hasta que haya esquema clínico local. |
| Advertir fuera de dominio calibrado | Implementado | `RegularGridInterpolator(bounds_error=True)` bloquea extrapolación. |

## Validación real

Las notas que requieren datos clínicos reales, calibraciones experimentales,
tablas publicadas digitalizadas o Monte Carlo de transporte se dejan como líneas
futuras. No se han inventado valores para `BSF`, transmisión de mesa ni
`f_skin`.
