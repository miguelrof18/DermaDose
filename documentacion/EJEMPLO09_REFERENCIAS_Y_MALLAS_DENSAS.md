# Mesa del ejemplo09 y mallas densas

El error «sistema de mesa isocenter diferente del configurado» aparece cuando el RDSR contiene posiciones relativas al isocentro y el perfil elegido describe posiciones legacy con un origen de mesa configurado.

## Adaptación

normalize conserva primero el triplete nativo del sistema configurado cuando está disponible en el RDSR. Si solo dispone de isocenter y el destino es legacy, calcula:

- posición en sala = isocenter_mm + matriz_origen × posición_RDSR;
- posición legacy = table_origin_mm + inversa(table_matrix) × posición_en_sala.

matriz_origen utiliza isocenter_table_matrix si se declara y, si no, los ejes elegidos mediante table_matrix. No se deduce una orientación de otro fabricante a partir de su nombre. La nueva matriz debe ser ortonormal. Esta adaptación expresa las coordenadas en los ejes seleccionados; no caracteriza por sí sola una instalación ni registra un paciente.

Se registran las coordenadas de origen, destino y sala en seleccion_entrada.json y en los eventos del informe. No se alteran el RDSR, el STL o el JSON originales. La validación se realiza antes de cargar la malla.

## Rendimiento

La malla utilizada en la comprobación contiene 2 042 996 triángulos. Una reproducción de cinco eventos requiere más de un millón de consultas de visibilidad. Se mantiene la resolución numérica de la malla:

- El motor RTree consulta cajas de segmentos de rayo en vez de una única caja grande. Las intersecciones finales se siguen resolviendo con el mismo núcleo de triángulos.
- Para mallas de al menos 200 000 caras, si existe PythonSlicer se puede utilizar vtkStaticCellLocator en un proceso local oculto. Se transfieren solo malla y rayos; no se abre una escena de Slicer. El índice se invalida al cambiar vértices o caras. Si no está disponible, se conserva RTree.
- DERMADOSE_PYTHONSLICER permite indicar la ruta del ejecutable. En Windows se busca también una instalación en AppData/Local/slicer.org. No se instala ningún paquete.
- Los resultados registran el motor y su versión. Las parejas candidatas no instrumentadas por VTK figuran como no disponibles, no como una medida de cero.
- Se liberan cachés grandes de topología y el trabajador nativo antes de los informes.
- La media local de 1 cm² es una métrica adicional. Por defecto no se calcula para mallas con más de 200 000 caras: figura null con estado not_computed_mesh_budget. La PSD puntual y los mapas se calculan sobre todas las caras. max_neighbourhood_faces permite cambiar el límite; 0 solicita la búsqueda completa sin límite.
- La imagen PNG 3D de mallas grandes es una vista de puntos muestreada e incluye el máximo. El mapa numérico conserva todas las caras y el HTML mantiene la malla completa. Ese HTML puede ser grande.
- Se muestra el avance por evento, sin añadir mensajes ni HTML de orientación del paciente.

## Comprobaciones y límites

Pasaron 67 pruebas de regresión antes de los últimos ajustes de recursos. Además se comparó el motor VTK con RTree en fuentes axiales, oblicuas e interiores; los identificadores de las primeras intersecciones coincidieron, incluyendo los puntos finales sobre la superficie.

La reproducción local de cinco eventos, perfil MRC200 y STL completo superó todos los eventos y llegó a las métricas. La conexión al equipo se interrumpió antes de confirmar la exportación final de ese caso grande. El límite de la métrica opcional y la vista PNG acotada se añaden para evitar esos costes; su comprobación automática se deja en la suite de pruebas.

Al revisar los archivos actuales, ejemplo09 contenía 0_anon.dcm (101 eventos), no intervencion_sintetica.dcm del mensaje original. Ese archivo tiene SID/SOD/SRP iguales a cero en el evento 91. Sigue rechazándose: no se han inventado distancias ni descartado su kerma.
