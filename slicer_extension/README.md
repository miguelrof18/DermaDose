# DermaDose Surface para 3D Slicer

Módulo Python independiente para extraer del TC un **contorno externo** como STL.
No detecta una capa histológica de piel y no ejecuta el cálculo dosimétrico.

## Instalación local

1. Instala una versión reciente de [3D Slicer](https://www.slicer.org/).
2. Arrastra `DermaDoseSurface/DermaDoseSurface.py` a la ventana de Slicer.
3. Elige **Add Python scripted modules to the application** y carga el módulo.
4. Abre **DermaDose → DermaDose Surface**.

La documentación de Slicer describe también la ruta **Edit → Application settings → Modules → Additional module paths**. No hace falta compilar ni publicar una extensión.

## Uso

1. Importa una serie TC DICOM a Slicer y selecciónala en el módulo.
2. Ajusta el umbral inicial (valor sugerido: -350 HU) y pulsa **Extraer superficie**.
3. Comprueba cortes axial, coronal y sagital y la vista 3D. Usa **Abrir Segment Editor para corregir** si la mesa, cables, inmovilizadores o fragmentos ajenos se han incluido; por ejemplo con Scissors. La mesa puede estar conectada al cuerpo y sobrevivir a la selección del componente mayor.
4. Pulsa **Revisar malla**, inspecciona sus avisos y marca la casilla tras la revisión visual.
5. Exporta. Se escriben `nombre.stl` y `nombre.json` (SHA-256, RAS/mm, matriz del TC, umbral y resultados de comprobación).

El STL binario lleva `SPACE=RAS` en su cabecera. Slicer trabaja internamente en RAS; los vértices se guardan en **milímetros**. El JSON no almacena identificadores del paciente, pero el STL conserva la anatomía: trátalo como dato sensible.

## Uso con DermaDose

El STL se puede seleccionar mediante `patient.mesh_path` en la configuración del equipo. **Antes de ejecutar el cálculo**, confirma que DermaDose interpreta el STL en RAS y milímetros y que la transformación entre el TC y la geometría de irradiación RDSR está determinada. Esta extensión no registra automáticamente TC y RDSR. Un contorno recortado por el campo del TC deja piel sin representación; el informe indica qué caras del volumen se tocan.

El módulo sustituye su segmentación anterior al volver a extraer. Guarda la escena de Slicer antes si deseas conservar una edición manual. La malla final depende de la resolución del TC, el umbral y las correcciones realizadas.

## Prueba en Slicer

En **Python console**:

```python
import unittest
from DermaDoseSurface import DermaDoseSurfaceTest
unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(DermaDoseSurfaceTest))
```

La prueba crea un volumen sintético, elimina una isla, rellena una cavidad, genera la superficie y verifica el STL y JSON. Debe ejecutarse dentro de Slicer; el `pytest` habitual de DermaDose no incluye las bibliotecas de Slicer.
