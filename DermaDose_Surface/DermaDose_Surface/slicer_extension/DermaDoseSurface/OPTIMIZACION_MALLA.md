# Optimización de mallas dentro de DermaDose Surface 1.6

El módulo busca la menor malla aceptable entre las alternativas ensayadas,
dentro de 100 000–200 000 caras. 150 000 es una referencia, no un requisito exacto.
No se afirma un óptimo global ni precisión dosimétrica a partir del número de caras.

## Implementado

1. **QEM con normales y conservación de volumen.** El error cuádrico evalúa el coste
   de colapsar aristas. La contribución de las normales penaliza cambios de forma;
   permite una distribución no uniforme según la geometría sin imponer una cuadrícula.
   No es el algoritmo de campo de curvatura de CGAL.
2. **Regularización de cuádricas probabilísticas de VTK**, cuando existe en la versión
   instalada (desde VTK 9.3). El parámetro es como máximo 0,05 y se reduce con la tolerancia
   geométrica. No aplica un filtro de suavizado a la anatomía.
3. **Alternativa con conservación de topología**, contrastada con QEM. Se ensayan
   presupuestos de 150k y 100k; 200k se ensaya cuando no hay una opción más pequeña válida.
   Para el par de comparación se utilizan los márgenes propios de cada malla.
4. **Selección por coste y forma.** Se selecciona el menor grupo de 500 caras que cumple
   los controles; entre opciones del mismo grupo se prefieren menos caras muy alargadas
   y menor desviación. Las alternativas siempre parten del mismo original.
5. **Calidad de triángulos**, usando la razón media
   `q = 4 sqrt(3) A / (a²+b²+c²)`: 1 para equiláteros, 0 para degenerados.
   Se guardan percentiles, fracción con q<0,1, área y longitud de arista.
   El percentil 1 debe ser como mínimo min(0,05, percentil 1 original).
6. **Control geométrico**: un componente, superficie cerrada, ausencia de degenerados,
   aristas no manifold y conservación de la característica de Euler. Cambio de volumen
   ≤1%, área ≤5%, desviación muestreada ≤tolerancia seleccionada en Slicer (1 mm por defecto).
   La desviación comprueba TODOS los vértices y TODOS los centroides de caras en ambos
   sentidos. No es una cota continua de Hausdorff ni detecta todas las autointersecciones.
7. **Relectura del STL**, comprobando integridad, número de caras y cambio adicional
   por serialización ≤0,001 mm. No se sobrescriben archivos existentes.
8. **Interfaz Slicer** con tolerancia ajustable, resumen de reducción, desviación,
   cambio de volumen y caras alargadas. Tanto extracción como nueva optimización tras
   editar la segmentación se ejecutan en procesos cancelables. Exportar una revisión
   obsoleta solicita revisar; no inicia una reducción larga en el hilo de interfaz.

Las mallas originales por debajo del margen se conservan sin subdividirlas.
Si ninguna opción cumple las condiciones, queda el original para revisión y se explica
el fallo. No se rellenan agujeros para satisfacer el número. La extracción previa del
contorno externo mantiene su tratamiento de cavidades internas; debe revisarse en los cortes.

## Métodos investigados y alcance

| Método | Decisión |
| --- | --- |
| QEM, error de normales, volumen y cuádricas probabilísticas | Integrado con VTK, ya disponible en Slicer; cada salida se valida |
| Remallado isotrópico adaptativo de CGAL | Investigado; no integrado. Combina colapsos, divisiones, giros, suavizado y reproyección; requiere otro backend y validación |
| Hausdorff acotada de CGAL | Investigada; no integrada. La comprobación VTK es muestreada y se identifica así |
| Reparación automática, Poisson, cierre de agujeros, suavizado global | No utilizados para forzar el margen: pueden cambiar la anatomía o topología |
| Refinamiento guiado por dosis | Requiere RDSR/registro y una validación aparte; no se presupone en un módulo de segmentación |

## Referencias primarias

- Garland y Heckbert, Surface Simplification Using Quadric Error Metrics (1997): https://www.cs.cmu.edu/~garland/quadrics/
- VTK vtkQuadricDecimation: https://vtk.org/doc/nightly/html/classvtkQuadricDecimation.html
- Regularización añadida en VTK 9.3: https://docs.vtk.org/en/v9.7.0/release_details/9.3.html
- CGAL, remallado con campos uniformes y adaptativos: https://doc.cgal.org/6.2/PMP_Remeshing/index.html
- CGAL, distancias de Hausdorff acotadas: https://doc.cgal.org/latest/Polygon_mesh_processing/

## Verificación y límites

Pruebas nativas VTK: reducción, búsqueda de densidad, forma, tolerancias, relectura STL,
conservación de un asa (toro), rechazo de superficies abiertas, exportación del par,
revisión mediante proceso separado sin modificar el original, y un caso de caras no
planas donde los vértices coinciden pero sus interiores se separan.

En una esfera sintética de radio 100 mm, 358 200 caras se redujeron a 100 000 con una
desviación muestreada aproximada de 0,012 mm. Es evidencia de funcionamiento, no de
exactitud para anatomías reales. La interfaz Qt/MRML sigue pendiente de ejecución en
Slicer; las pruebas nativas no la sustituyen. Mantener la PR como borrador hasta verificarla.

El número óptimo para DOSIMETRÍA debe justificarse comparando PSD en las dos mallas
independientes de 150k/300k con el mismo RDSR, configuración y registro.
