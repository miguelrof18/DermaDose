# Ejecutar ejemplos con XCAT o malla del paciente

En cada carpeta `ejemplos rdsr/ejemploXX/`, hermana del programa, guarda el RDSR .dcm y el JSON del equipo. Para usar una superficie de TC, añade su .stl con cualquier nombre. Para XCAT no hace falta copiar una malla: se usa la incluida en el programa.

Al ejecutar EJECUTAR_EJEMPLO.bat:

1. Introduce el número de ejemplo.
2. Elige **1: XCAT** o **2: TC del paciente (STL local)**.
3. Si hay varios archivos del tipo elegido, selecciona uno.
4. Consulta el resumen de colocación; el cálculo continúa automáticamente.

Los JSON de metadatos de Slicer no se confunden con equipos. No se requiere paciente.json, ni se utiliza aunque exista; tampoco se utilizan por ahora los metadatos acompañantes del STL. La colocación procede del bloque patient del equipo. Los campos ausentes se completan con mm, orientación identidad, traslaciones nulas, resolución de 5 mm y espesor de 220 mm. No se centra automáticamente la superficie. El perfil y las tablas espectrales del equipo se conservan.

## Resumen inicial legible

En consola se muestran las direcciones de los ejes del STL, las traslaciones en centímetros y las inclinaciones adicionales. La ficha resultados/colocacion_malla.html presenta la misma información con tarjetas de colores y una tabla; la matriz se puede desplegar como detalle técnico. Con el .bat se abre antes de iniciar el cálculo y el mapa 3D se abre al terminar.

La convención del motor es +X izquierda, +Y cabeza, +Z anterior. Por ejemplo, con la matriz XCAT [[1,0,0],[0,0,1],[0,-1,0]], el eje +X del STL apunta a izquierda, +Y a posterior y +Z a cabeza. Una traslación [0,-500,40] se describe como **50 cm hacia pies y 4 cm hacia anterior**.

Estas etiquetas explican la transformación configurada: no identifican la anatomía de una malla arbitraria ni confirman su registro real. La traslación de la malla se aplica antes de las inclinaciones adicionales; el desplazamiento adicional respecto a la mesa y sus movimientos del RDSR se aplican después. La colocación se registra como provisional sin bloquear la ejecución. El resumen no deduce posiciones clínicas como cabeza primero o decúbito a partir del STL.

## Comandos

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier.example_runner "..\ejemplos rdsr\ejemplo07" --interactive
.\.venv\Scripts\python.exe -m skindose_verifier.example_runner "..\ejemplos rdsr\ejemplo07" --mesh xcat
.\.venv\Scripts\python.exe -m skindose_verifier.example_runner "..\ejemplos rdsr\ejemplo07" --mesh superficie.stl
```

--prepare-only valida y guarda la configuración y el resumen sin calcular dosis. --open-viewer abre la ficha inicial y el mapa final. Los nombres explícitos se interpretan desde la carpeta del ejemplo. El comando clásico skindose_verifier run mantiene su comportamiento.

## Tablas y resultados

Las referencias relativas conocidas a skindose_verifier/data/derived_mrc200_nominal, derived_azurion_spekpy y derived_azurion_qc_conditional se resuelven contra los datos del programa si el JSON se ha copiado desde config/. Una ruta personalizada inválida produce error y no cambia de tablas.

Los originales no se reescriben. resultados/configuracion_efectiva.json conserva los datos usados y resultados/seleccion_entrada.json registra rutas, hashes, selección XCAT/paciente y colocación provisional (registration_status: not_established). Al repetir se actualizan las salidas de esa carpeta.

## Distancias nulas o ausentes en un evento Azurion

El modo habitual rechaza distancias inválidas e identifica el índice del evento, su UID, sus distancias y el kerma asociado. Un evento con dosis y geometría inválida no se elimina del cálculo.

Para ensayos con un perfil Azurion que ya permita hipótesis de investigación, se puede solicitar expresamente:

```powershell
.\EJECUTAR_EJEMPLO.bat --research-distances
```

También está disponible --research-distances al invocar directamente skindose_verifier.example_runner. Solo sustituye distancias ausentes, no finitas o no positivas por los valores de geometry del JSON correspondiente al plano. Conserva todas las distancias válidas y el kerma de cada evento. No interpola ni copia geometría de eventos vecinos.

La sustitución no recupera una medida real. La opción solo afecta a la configuración efectiva de esa ejecución; no modifica el RDSR ni el JSON de entrada. La consola y seleccion_entrada.json muestran cada evento afectado y, para cada distancia, el valor original, el usado y su fuente. El informe por evento conserva la marca de sustitución de investigación.

La comprobación geométrica se realiza antes del cálculo sobre la malla. Superarla no valida por sí sola la colocación, los factores dosimétricos ni el resultado final.
