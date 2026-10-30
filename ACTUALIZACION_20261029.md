# Actualización integrada — 29 de octubre de 2026

Se integraron las correcciones de la revisión en esta carpeta concreta. Se conservaron los archivos originales, mallas, configuraciones, resultados y entornos virtuales incluidos en el ZIP.

## Verificación sobre esta copia

- 119 pruebas del programa superadas, sin fallos.
- 6 pruebas independientes VTK del módulo Surface superadas. El archivo del módulo coincide con el revisado y no necesita reemplazo.
- Las mallas originales aportadas resuelven los siete bloqueos documentados en la revisión anterior. Los resultados antiguos de la memoria describen aquella copia; el resultado de esta integración está en los registros adjuntos.
- No se ejecutó la interfaz completa en Windows ni una simulación clínica del Azurion.

## Cambios

Se incorporaron los controles de filtración, el rechazo por defecto de sustituciones sintéticas, la política de interpolación combinada, comprobaciones del perfil Azurion en llamadas Python, búsqueda del máximo promediado, validación de incertidumbres y referencias, y empaquetado de tablas.

Se añadió la tabla espectral condicionada de 910 filas y su configuración alternativa `config/equipo_azurion_espectral_condicionado.json`. Su Al efectivo y ángulo son hipótesis condicionadas; no son una caracterización completa. No se ha activado esta configuración automáticamente ni modificado tus perfiles actuales.

La memoria PDF, fuentes LaTeX e informe previo están en `documentacion/revision_20261029/`, junto con los registros de esta integración.

## Cómo usarla

1. Descomprime la carpeta DermaDose completa. Conserva una copia de tu carpeta anterior.
2. Mantén `ejemplos rdsr` al mismo nivel que DermaDose: los ejecutables BAT buscan los ejemplos en esa carpeta hermana, que no estaba dentro del ZIP recibido.
3. Usa tus ejecutables habituales y selecciona la configuración correspondiente. Los resultados guardados anteriores se conservan, pero no se recalcularon.
4. Si necesitas recrear el entorno virtual porque ha cambiado su ruta, desde la carpeta DermaDose usa Python compatible y ejecuta:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[test,spectral]"
.\.venv\Scripts\python.exe -m pytest tests
```

Los entornos Windows incluidos se han conservado; las pruebas se realizaron con un entorno Linux independiente, no con esos ejecutables Windows.
