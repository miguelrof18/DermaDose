# Ejecutar ejemplos con XCAT o malla del paciente

En cada carpeta de ejemplo guarda el RDSR .dcm y el JSON del equipo. Para usar una superficie de TC añade su .stl con cualquier nombre. XCAT utiliza la malla incluida en el programa.

Al ejecutar EJECUTAR_EJEMPLO.bat:

1. Introduce el número de ejemplo.
2. Elige **1: XCAT** o **2: TC del paciente (STL local)**.
3. Si hay varios archivos del tipo elegido, selecciona uno.
4. Se inicia el cálculo y se abre el mapa 3D al terminar.

No se requiere paciente.json. La colocación del bloque patient del equipo se sigue aplicando internamente sin mostrar mensajes ni generar o abrir una ficha de orientación. Al ejecutar de nuevo el ejemplo se elimina el antiguo resultados/colocacion_malla.html si existe. No se deduce automáticamente la colocación real de una superficie a partir de su STL.

El JSON de metadatos de Slicer no se confunde con el JSON del equipo. El RDSR, el JSON y el STL originales no se reescriben.

## Distancias del RDSR

No existe la opción --research-distances. Los eventos Azurion con distancias ausentes, no finitas o no positivas detienen el cálculo con un mensaje que identifica el evento y sus datos inválidos. No se sustituyen por valores del JSON ni se omite su dosis. Una antigua configuración con invalid_distance_policy: configured_research tampoco activa sustituciones.

Esta retirada afecta a la sustitución de distancias; conserva las hipótesis dosimétricas declaradas en los perfiles de equipo existentes.

## Comandos

```powershell
.\EJECUTAR_EJEMPLO.bat
.\.venv\Scripts\python.exe -m skindose_verifier.example_runner "..\ejemplos rdsr\ejemplo07" --mesh xcat
.\.venv\Scripts\python.exe -m skindose_verifier.example_runner "..\ejemplos rdsr\ejemplo07" --mesh superficie.stl
```

--prepare-only comprueba y guarda la configuración sin calcular dosis. --open-viewer abre únicamente el mapa de dosis final.

## Tablas y resultados

Se mantienen las rutas conocidas de tablas espectrales y los factores del equipo elegido. resultados/configuracion_efectiva.json y resultados/seleccion_entrada.json conservan la configuración efectiva y la selección de archivos. Las coordenadas necesarias para reproducir el cálculo permanecen en esos datos técnicos. Al repetir se actualizan las salidas de esa carpeta.
