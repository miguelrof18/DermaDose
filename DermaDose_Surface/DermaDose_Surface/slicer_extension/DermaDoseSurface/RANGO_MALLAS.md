# Exportación de superficies: 150k y comparación 300k

DermaDose Surface 1.5 busca 150 000 triángulos y acepta 100 000–200 000.
Si el original es más pequeño, se conserva sin subdividirlo.
No se usa suavizado ni relleno de agujeros para alcanzar el número.
El relleno de cavidades internas pertenece a la extracción previa del contorno externo del TC,
no a la simplificación; revise siempre el contorno, la mesa y los límites de imagen.

La reducción exige superficie cerrada, un componente, ausencia de aristas no manifold
y caras degeneradas, y conservación de la característica de Euler.
Se comprueba además una distancia bidireccional entre vértices y la superficie opuesta
de como máximo 1 mm y un cambio de volumen de como máximo 1%.
Estas tolerancias son controles geométricos iniciales, no una validación de precisión dosimétrica.
La distancia entre vértices y superficie no es una cota continua de Hausdorff.

Si ninguna reducción cumple las condiciones se conserva el original para revisión,
se explica el fallo y no se exporta silenciosamente una malla de millones de caras.
El STL se vuelve a leer y validar antes de publicarlo. No se sobrescriben archivos.

## Comparación

El botón «Exportar comparación 150 000 / 300 000» genera ambas mallas
independientemente desde el mismo original, en segundo plano y con cancelación.
Para 300k el margen es 270 000–330 000. Requiere un original de al menos 300 000 caras;
no interpola triángulos para simular una referencia más detallada.
Se escriben dos STL, sus JSON y comparison.json con el hash del original.
El proceso conserva coordenadas RAS/mm. Editar la segmentación invalida la caché reducida.

Con una configuración de malla cuyo registro ya se haya revisado:

```powershell
python -m skindose_verifier.mesh_convergence informe.dcm --config equipo.json --low surface_150k.stl --high surface_300k.stl --out comparison_psd.json
```

El cálculo usa el mismo RDSR, geometría, factores y registro. El informe incluye ambas PSD puntuales,
diferencia en Gy, porcentaje respecto a 300k y desplazamiento del máximo.
No modifica distancias faltantes ni introduce factores sustitutivos.
La comparación debe repetirse con casos representativos; ni 150k ni 300k garantizan por sí solos la precisión.

## Estado de verificación

Pruebas automáticas de geometría VTK: reducción, topología, tolerancias,
relectura STL, rechazo de superficies abiertas y exportación del par.
La validación de interfaz de Slicer y la comparación PSD con TC/RDSR reales requieren ejecución
en el entorno del proyecto; el entorno local dejó de estar disponible durante este cambio.
