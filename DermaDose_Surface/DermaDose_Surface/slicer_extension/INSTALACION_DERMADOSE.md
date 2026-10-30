# Instalación en tu proyecto DermaDose

1. Copia la carpeta slicer_extension de este paquete en la raíz de tu repositorio DermaDose, junto al README del proyecto. No la copies dentro del entorno virtual.
2. Abre 3D Slicer. En Edit > Application settings > Modules > Additional module paths añade la carpeta DermaDose/slicer_extension/DermaDoseSurface (la que contiene DermaDoseSurface.py). Reinicia Slicer y busca DermaDose Surface en el selector de módulos.
3. Carga la serie TC usando el módulo DICOM. Abre DermaDose Surface y selecciona el volumen.
4. Pulsa Extraer superficie. El umbral inicial -350 HU es orientativo. Revisa los cortes y la vista 3D; corrige mesa y objetos ajenos en Segment Editor. Vuelve al módulo, pulsa Revisar malla, revisa los avisos y marca la casilla.
5. Exporta como paciente.stl. Se crea también paciente.json. Guarda ambos, por ejemplo, en una carpeta de datos local de DermaDose. No subas anatomía de pacientes a GitHub.

## Conexión con el cálculo de DermaDose

Este módulo se ejecuta dentro de Slicer, no con python ni pip en el entorno de DermaDose. Su salida es la malla; no instala automáticamente un nuevo lector en el motor de dosis.

La guía original del módulo indica patient.mesh_path como opción para seleccionar el STL. No se ha podido comprobar aquí el esquema actual del repositorio: confirma ese nombre en la configuración y en el cargador de paciente de tu versión antes de editarlo. Si existe, apunta esa opción al STL exportado. Conserva el JSON junto al STL; no se ha comprobado que el motor lo lea automáticamente.

El STL se exporta en RAS y milímetros, con el origen del TC. No basta con seleccionar el archivo: el motor debe admitir STL y usar las mismas unidades y coordenadas. La matriz ijk_to_ras_mm del JSON describe el TC, NO el registro entre TC y equipo. Es necesario determinar esa transformación antes de calcular dosis. Si el motor usa LPS, RAS a LPS cambia el signo de X e Y; eso por sí solo no registra al paciente respecto al equipo. Si usa metros, también hay que convertir mm a m.

XCAT50 puede seguir como configuración predeterminada. Usa una configuración separada para la malla TC una vez resuelta la transformación; copiar esta carpeta no cambia el fantoma actual.

## Estado

Sintaxis Python comprobada. Slicer no está instalado en este entorno: la interfaz, la prueba sintética y la exportación deben probarse dentro de Slicer. Consulta README.md para ejecutar la prueba incluida. La superficie es el contorno externo del volumen disponible, no una capa histológica de piel.
