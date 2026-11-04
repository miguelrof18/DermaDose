# DermaDose CT independiente — 2.0.0

Nuevo flujo de segmentación y exportación para ejecutar DENTRO de 3D Slicer.
Es independiente del módulo DermaDose Surface instalado: no hay que reemplazar módulos.
Mantenga dermadose_ct.py y surface_processing.py en la misma carpeta.
No instale el paquete PyPI llamado slicer: no es la aplicación 3D Slicer.

## Abrir sin escribir rutas

1. Descomprima el ZIP en una carpeta nueva.
2. Abra 3D Slicer y cargue el TC (volumen escalar, valores HU).
3. Abra la consola Python de Slicer (View → Python Interactor).
4. Ejecute estas líneas. El selector permite localizar el archivo recién descomprimido:

```python
import runpy, qt
ruta = qt.QFileDialog.getOpenFileName(None, 'Seleccionar dermadose_ct.py', '', 'Python (*.py)')
dd = runpy.run_path(ruta, run_name='__main__') if ruta else None
```

En la ventana «DermaDose CT independiente 2.0.0-standalone»:

1. Seleccione el TC. Umbral inicial: −350 HU. Es orientativo y ajustable.
2. Pulse «Segmentar cuerpo». Aparece el segmento en los cortes, sin generar aún una malla grande.
3. Pulse «Revisar y corregir en Segment Editor». Compruebe mesa, ropa, brazos y cobertura.
   Use Scissors para retirar la mesa si está conectada. Revise que los brazos no hayan
   desaparecido al seleccionar el componente más grande. Evite pulsar Show 3D sobre
   la segmentación de alta resolución: la vista 3D reducida se prepara en el paso siguiente.
4. Pulse «Preparar malla 3D». El cálculo se realiza en PythonSlicer, con progreso y cancelación.
5. Tras revisar la vista 3D pulse «Exportar STL + JSON + segmentación». Elija un nombre nuevo.

Si cambia el segmento o el TC, vuelva a preparar la malla antes de exportar.
Si la reducción falla, la segmentación sigue en la escena: puede corregirla y guardarla
con el botón Save de Slicer. No se destruye porque no se alcance el recuento de caras.

## Qué hace

- Umbral de tejido, componente 3D principal (6 vecinos) y relleno de cavidades internas.
- Opcionalmente rellena huecos por corte de la rejilla de adquisición. Esto ayuda a
  obtener el contorno externo incluso si la vía aérea conecta con el exterior. Puede
  rellenar huecos que se desee conservar: está identificado y puede desactivarse.
- No usa apertura, cierre ni suavizado global de la anatomía por defecto.
- Construye la superficie solamente en la caja que contiene el segmento, conservando
  su desplazamiento IJK y la matriz IJK→RAS completa, también en TC oblicuos.
- Reduce hacia 150 000 caras, aceptando 100 000–200 000. Las mallas menores se conservan
  sin añadir caras. Usa QEM y una alternativa topológica con validación de desviación,
  volumen, área, topología y calidad de triángulos. El código está en surface_processing.py.
- Exporta STL binario RAS/mm, JSON de auditoría y .seg.nrrd editable. Relee el STL y evita
  sobrescribir archivos. Esto no establece el registro con el equipo intervencionista.
- La reducción se intenta una vez y solo prueba alternativas cuando hace falta.
  La preparación, importación MRML, copia de datos y guardado todavía pueden producir
  pausas breves; no se promete que todos los TC sean rápidos ni que todos pasen los controles.

## Limitaciones que hay que revisar

Un TC no incluye piel fuera del campo adquirido. La superficie se cierra en sus límites;
esas tapas artificiales no son piel real. El programa las señala en la auditoría.
El componente mayor no distingue paciente y mesa si están conectados. Tampoco conserva
por sí solo brazos separados del tronco. Recorte con Crop Volume o corrija en Segment Editor.
El perfil limita la entrada a 160 millones de vóxeles para acotar memoria; no es un umbral
universal de hardware. Es un contorno corporal externo, no una capa histológica de piel.
150 000 caras no garantizan precisión dosimétrica; compare PSD en mallas de distinta resolución.

## Extensiones y procedimiento documentado

La receta de Andras Lasso para superficies de piel recomienda Threshold → Islands /
Keep largest island → Wrap Solidify / Outer surface. Está documentada para CT/MRI;
el ejemplo de intensidad de MRI no se debe copiar como umbral HU de un TC.

- Receta: https://lassoan.github.io/SlicerSegmentationRecipes/SkinSurface2/
- Extensión SurfaceWrapSolidify: https://github.com/sebastianandress/Slicer-SurfaceWrapSolidify
- Remove CT Table, SlicerSandbox: https://github.com/PerkLab/SlicerSandbox/blob/master/RemoveCtTable/RemoveCtTable.py
- Efectos y parámetros de Segment Editor: https://slicer.readthedocs.io/en/latest/developer_guide/modules/segmenteditor.html
- Segmentación desde Python: https://github.com/Slicer/Slicer/blob/main/Docs/developer_guide/script_repository/segmentations.md
- Artículo sobre reconstrucción de superficie de cabeza desde TC (2025):
  https://pmc.ncbi.nlm.nih.gov/articles/PMC12109114/
- Pinter et al., «Polymorph segmentation representation for medical image computing»,
  Computer Methods and Programs in Biomedicine 171 (2019), 19–26:
  referencia enlazada por https://slicer.readthedocs.io/en/5.8/user_guide/modules/segmenteditor.html

El script nuevo implementa la vía de umbral/componente/relleno sin instalar extensiones.
El relleno por vóxeles NO es el algoritmo Wrap Solidify y no se afirma equivalencia.
Si hay conexiones difíciles o cavidades complejas, instale SurfaceWrapSolidify desde
Extensions Manager, reinicie Slicer y aplique su efecto Outer surface en Segment Editor
sobre el segmento. Después vuelva a preparar la malla con este script.
Remove CT Table también necesita Wrap Solidify según su implementación; no resuelve
universalmente el contacto directo paciente-mesa.

## Validación y autoprueba en Slicer real

Se han superado siete pruebas locales con numpy/scipy/VTK: selección corporal, cavidades,
mesa separada y conectada, rejilla oblicua con recorte, tapas por truncación, entradas inválidas
y ejecución del trabajador que produce STL y metadatos. El Python se compila correctamente.
No se ha podido ejecutar esta ventana en una instalación completa de Slicer en esta sesión.

Tras abrir el script, ejecute en la consola de Slicer:

```python
dd['self_test']()
```

Crea un TC sintético pequeño y prueba segmentación, QProcess, vista 3D y los tres archivos
exportados, conservando los nodos previos de su escena. Elimina solo los nodos de prueba.
Si termina correctamente imprime DERMADOSE_SELF_TEST_OK. Si falla, copie el traceback
completo y el contenido del cuadro de estado. La autoprueba está incluida pero aún no ha
sido ejecutada en Slicer real por el autor de esta entrega.
