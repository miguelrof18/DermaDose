# Revisión científica y técnica de DermaDose — 29 de octubre de 2026

Se corrigió la copia de código disponible y se actualizó la memoria LaTeX. Esta entrega no acredita validación clínica ni una auditoría exhaustiva cerrada: hay siete pruebas bloqueadas por recursos ausentes y no se ha ejecutado la aplicación completa en Windows/3D Slicer.

## Hallazgos corregidos

| Hallazgo | Corrección | Evidencia |
|---|---|---|
| Sustitución automática de geometría ante equipo distinto | `allow_synthetic_device_override` booleano explícito, falso por defecto; prohibido con perfil Azurion | Regresión de rechazo por defecto |
| Sustitución automática de factores publicados fuera de dominio | `allow_historical_factor_fallback` explícito; sustitución global, nunca mezcla por evento; prohibido en Azurion | Regresión fuera de dominio |
| Media injustificada de espesores mínimo/máximo | Rechazo si difieren; tolerancia numérica explícita | Caso de filtro no uniforme |
| Etiqueta sintética permitía ignorar Al añadido | La política `reject_additional` se aplica siempre | Regresión de elusión |
| Generador sin política de interpolación conjunta | Dataset nuevo versión 1.1 declara `combined_C`; conserva C/B como factor efectivo | Generación real con SpekPy |
| API Python podía omitir comprobaciones del perfil Azurion | Comprobaciones en `calculate`; configuración local copiada | Pruebas de rechazo y no mutación |
| Máximo promediado podía excluir centro con dosis nula | Búsqueda sobre todos los centros | Caso discreto de media 6,6 Gy |
| Promedio en elipse sin resolución longitudinal | Devuelve no disponible | Inspección de condición |
| Incertidumbres RDSR negativas/no finitas | Validación al construir configuración | Tres regresiones |
| Referencias comerciales NaN/infinito aceptadas | Rechazo explícito | Dos regresiones |
| Incertidumbre cuadrática podía confundirse con incertidumbre de PSD | Metadato declara máximos de componentes y covarianza nula | Inspección de salida |
| Dataset ampliado ausente de distribución instalable | Inclusión de CSV/JSON en package-data | Wheel construido e inspeccionado |
| Procedencia anatómica atribuida genéricamente | Descripción neutral en el motor | Inspección |

Se preservaron las tablas existentes. No se han modificado resultados históricos para aparentar una mejora de exactitud. Los datasets antiguos conservan su política original; regenerar una tabla con el constructor nuevo cambia su versión y política de interpolación interior.

## Verificación

- Batería: **119 passed, 7 failed**, 126 casos. Los siete fallos permanecen visibles; necesitan `xcat_wholebody_body.stl`, `pyskindose_adult_female.stl` o configuración original asociada. No son siete errores dosimétricos demostrados ni siete casos aprobados.
- Regresiones nuevas: 13 casos aprobados; contra la copia anterior, 12 fallos y un caso aprobado. Se conservaron los registros.
- Se ajustaron pruebas antiguas que esperaban sustituciones implícitas: ahora deben autorizarlas. El caso de filtro homogéneo usa extremos iguales, y se añadió rechazo específico del no homogéneo.
- Reproducción espectral: 78 + 182 espectros, 390 + 910 filas; segunda implementación de interpolación de productos con SciPy y suma explícita. Diferencias relativas máximas de C: 9,992×10⁻¹⁶ y 8,882×10⁻¹⁶. HVL idénticas numéricamente. Hashes de CSV correctos.
- Se usa el mismo generador SpekPy y los mismos datos monoenergéticos: esto comprueba álgebra y reproducibilidad, no constituye una validación experimental independiente.
- Wheel construido; presencia comprobada de las dos tablas derivadas.
- Memoria: 106 páginas y 32 figuras vectoriales TikZ/PGFPlots; compilación sin referencias indefinidas ni cajas desbordadas; inspección visual de páginas modificadas.

## Correcciones de la memoria

Se precisó la conversión de kerma en aire libre a kerma superficial en agua y la aproximación necesaria para interpretarla como dosis. Se corrigieron filtración, cobertura histórica frente a ampliada, estado de observaciones HVL, significado de `thickness_mm`, transformaciones y colocación de malla. Se añadieron los booleanos de autorización y `sweep_samples` al anexo JSON y un capítulo de revisión con evidencia y límites.

La semántica de filtros se contrastó con DICOM PS3.17 UUUU.4.7: https://dicom.nema.org/medical/dicom/current/output/chtml/part17/sect_uuuu.4.7.html
El formalismo se contrastó con el informe aportado KS-ASF-201101-IR, sección 3.3, ecuación 13. La memoria conserva las referencias bibliográficas completas.

## Uso del código entregado

Descomprime el ZIP en una carpeta nueva. Contiene `dermadose/`, el módulo Slicer y pruebas, además de evidencia y un parche respecto a la copia inicial. No contiene TC ni datos clínicos ni mallas de fantoma que no estaban disponibles. No se ha instalado en tu ordenador ni enviado a GitHub.

En PowerShell, dentro de `dermadose`:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -e ".[test,spectral]"
.\.venv\Scripts\python.exe -m pytest tests/test_audit_regressions.py
.\.venv\Scripts\python.exe -m skindose_verifier run "RUTA_RDSR.dcm" --config "RUTA_CONFIG.json" --out "salida" --stage skin_dose
```

Comprueba la sintaxis disponible con `python -m skindose_verifier run --help` y conserva tus configuraciones/mallas originales. Para Azurion, usa el perfil explícito y rutas de dataset/malla existentes. No actives las banderas sintéticas para resolver carencias de caracterización del equipo real.

## Límites que permanecen

No se ha demostrado compatibilidad universal con los RDSR del fabricante, fidelidad anatómica de segmentación/simplificación, ausencia de autointersecciones, registro clínico, independencia de parámetros de incertidumbre o corrección del transporte disperso. Mesa escalar no incluye endurecimiento, y factores centrales de agua no resuelven oblicuidad y bordes en anatomía. La revisión de transcripción completa de todas las tablas bibliográficas y la campaña externa de validación siguen abiertas. El coste de evaluar todos los centros del promedio puede aumentar; no se ha medido rendimiento con una malla clínica de 150.000 caras.
