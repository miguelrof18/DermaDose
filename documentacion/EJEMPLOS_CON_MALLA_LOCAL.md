# Ejecutar ejemplos con la malla elegida

Cada carpeta `ejemplos rdsr/ejemploXX/`, hermana del programa, contiene:

- Un RDSR `.dcm`.
- Un JSON de equipo, por ejemplo `equipo_mrc200_conjunto_investigacion.json`.
- Una malla `.stl`, con cualquier nombre. Puede ser una copia del XCAT o una superficie segmentada.
- Opcionalmente el JSON de metadatos exportado por Slicer junto al STL.
- Opcionalmente `paciente.json`, con unidades y registro geom�trico propios del caso.

Ejecute `EJECUTAR_EJEMPLO.bat` e introduzca el n�mero. Un �nico archivo de cada tipo se selecciona autom�ticamente. Si hay varios RDSR, equipos o STL, el lanzador pregunta cu�l usar. Solo los JSON con version, manufacturer, model, geometry y factors se consideran equipos; el JSON de Slicer no se confunde con ellos.

El STL local sustituye la ruta de malla del JSON de equipo. Para usar XCAT puede copiar su STL a la carpeta, o invocar `--mesh xcat`. Si no existe ning�n STL local, el lanzador interactivo ofrece el XCAT incluido. No se cambia el modelo espectral ni se elige otro equipo de forma silenciosa.

## Colocaci�n del paciente

STL no contiene unidades ni registro con el arco. Los metadatos RAS de Slicer tampoco establecen d�nde estaba el paciente respecto al equipo. Si no existe paciente.json, se muestra el bloque patient del equipo y se pide aceptar expl�citamente sus unidades, orientaci�n y posici�n. Esas coordenadas pueden ser las del XCAT: no deben aceptarse a ciegas para una segmentaci�n de TC.

Para separar equipo y paciente, use `paciente.json` con este esquema, sustituyendo los valores por los del caso verificado:

```json
{
  "patient": {
    "mesh_unit": "mm",
    "mesh_axes_matrix": [[1,0,0],[0,1,0],[0,0,1]],
    "mesh_origin_mm": [0,0,0],
    "origin_on_table_mm": [0,0,0],
    "resolution_mm": 5,
    "mesh_max_edge_mm": null,
    "thickness_mm": 220
  }
}
```

Los ceros, la identidad y 220 mm son ilustrativos, no un registro de paciente. Las coordenadas se transforman como `p_local = escala * matriz * p_STL + mesh_origin_mm`; despu�s intervienen origin_on_table_mm y la geometr�a de mesa del motor. No se centra ni rota autom�ticamente una malla. mesh_max_edge_mm controla el refinamiento triangular; resolution_mm no fija el tama�o de tri�ngulo de un STL.

El JSON de equipo conserva los factores y la geometr�a del aparato. El bloque patient de paciente.json reemplaza el del equipo: incluya tambi�n el espesor si lo necesitan sus factores. El hash y unidades del metadato acompa�ante se comprueban cuando existen. La revisi�n de integridad de la malla contin�a en el motor habitual.

## Tablas espectrales

Al copiar un JSON desde config/ al ejemplo, las referencias relativas conocidas a `skindose_verifier/data/derived_mrc200_nominal`, `derived_azurion_spekpy` o `derived_azurion_qc_conditional` se resuelven contra los datos del programa si la ruta original ya no existe. Las rutas personalizadas inv�lidas producen error: no se sustituye una tabla por otra.

El perfil MRC200 conserva sus hip�tesis nominales de 2,5 mm Al equivalente y 11 grados. La selecci�n autom�tica del STL no caracteriza la instalaci�n ni altera su espectro.

## Comandos expl�citos

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier.example_runner "..\ejemplos rdsr\ejemplo07" --interactive
.\.venv\Scripts\python.exe -m skindose_verifier.example_runner "..\ejemplos rdsr\ejemplo07" --config equipo.json --mesh superficie.stl --patient paciente.json
.\.venv\Scripts\python.exe -m skindose_verifier.example_runner "..\ejemplos rdsr\ejemplo07" --mesh xcat --use-configured-placement
```

`--use-configured-placement` acepta expresamente la colocaci�n existente; `--prepare-only` prepara y valida la configuraci�n sin calcular dosis. Los nombres expl�citos se interpretan desde la carpeta del ejemplo. El comando cl�sico `skindose_verifier run` mantiene su comportamiento.

## Resultados y trazabilidad

Los originales nunca se reescriben. `resultados/configuracion_efectiva.json` registra las rutas resueltas y `resultados/seleccion_entrada.json` conserva la configuraci�n original, su hash, STL, hash, metadatos y origen de la colocaci�n. Los informes habituales quedan en resultados/. Al repetir se actualizan las salidas de esa carpeta; no representa un historial de ejecuciones.
