# Configuración del equipo

Usa un perfil por equipo y revisión. Mantén las transformaciones específicas del paciente fuera de los perfiles que se reutilizan en el lote.

| Bloque | Función |
|---|---|
| Identificación y versión | Fabricante, modelo y estado de caracterización. |
| `geometry` | Distancias de referencia, plano de colimación, signos y orden angular, origen y ejes de mesa, dimensiones del tablero. |
| `azurion_profile` | Reglas de lectura y recuperación de geometría del Azurion. La referencia de mesa desde el primer evento es una hipótesis provisional. |
| `patient` | Unidad, resolución y colocación de la superficie cuando corresponde. En la colocación automática el programa genera la transformación a partir del STL, el JSON de Slicer y el RDSR. |
| `factors` | Corrección del kerma y transmisión de mesa. Documenta la procedencia y la incertidumbre de las medidas. |
| `surface_factors`, `dose_medium` | Modelo de conversión a dosis, filtración, tablas y dominio de validez. |
| `tube_characterization`, `tube_output` | Validaciones del conjunto del tubo y rendimiento de referencia, cuando se configura. |

Verifica con un fantoma la lateralidad, los signos angulares y los desplazamientos de mesa antes de interpretar la colocación como representativa de la intervención.
