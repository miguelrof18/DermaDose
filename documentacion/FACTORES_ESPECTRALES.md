# Factores espectrales derivados para DermaDose

## Método

Módulo `skindose_verifier.spectral`. Implementa la ecuación 13 del informe KS-ASF-201101-IR: C(Q)=sum[K_i C(E_i)]/sum[K_i], con K_i proporcional a E_i*fluencia_diferencial_i*(mu_en/rho)_aire_i*ancho_bin_i y C(E)=B(E)*f(E). El producto se calcula primero en los nodos monoenergéticos publicados y luego se interpola dentro de su dominio. No se usa energía media ni se multiplican promedios independientes.

Para conservar la interfaz del motor se almacena B(Q)=<B(E)>_K y un factor efectivo f_eff=C(Q)/B(Q). Su producto reproduce exactamente C(Q) en cada nodo generado. f_eff NO es el cociente promediado por separado mediante la ecuación 12 ni un nuevo dato publicado por Benmakhlouf. La tabla se identifica como CALCULADA y no sustituye las tablas originales. Al interpolar entre nodos del dataset derivado, el proveedor existente interpola B y f_eff separadamente: su producto es una aproximación a C entre nodos. Refinar los ejes y estudiar convergencia si se requiere precisión adicional.

SpekPy 2.5.4 genera espectros de tungsteno. Se usan sus coeficientes de absorción de energía del aire NIST, no los de atenuación total. Se calcula HVL1 y HVL2 en Al y su cociente. Se registra versión, modelo físico, parámetros y hash del espectro y del CSV monoenergético original.

## Límites y comprobación

Tabla monoenergética disponible: 5–150 keV y campos cuadrados de lado 5–35 cm. No se extrapola energía ni campo. Los bins exteriores se aceptan solo si aportan como máximo 0,0001 del kerma total (0,01 %); se registra la fracción eliminada y la renormalización. El API permite tolerancia 0 para exigir cobertura íntegra. No se rellenan energías inferiores a 5 keV con valores inventados.

Antes de aplicar al equipo: contrastar el espectro con HVL1 medida y, cuando esté disponible, HVL2 y otras observaciones independientes. Una HVL compatible no identifica por sí sola el espectro ni valida los factores. El módulo no ajusta automáticamente parámetros para forzar acuerdo.

El benchmark incluido compara cuatro nodos clínicos del informe: diferencias en factor combinado alrededor de 0,06–0,14 % en la primera ejecución. El archivo de resultados registra también diferencias de HVL1 y B. Es una comprobación parcial de consistencia, no reproducción de todo el informe ni validación clínica. Los generadores de espectros pueden diferir y el producto de los factores separados de la tabla clínica no tiene por qué coincidir exactamente con la ecuación 13.

## Instalación y generación

Desde la raíz de DermaDose:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-spectral.txt
.\.venv\Scripts\python.exe -m skindose_verifier.spectral benchmark --anode-angle 12 --out comparacion_nueva.json
.\.venv\Scripts\python.exe -m skindose_verifier.spectral build --out tablas_espectrales_nuevas --kvps 46,50,60,70,80,90,100,110,120,125,130,140,150 --al 3.1,4.1 --cu 0,0.1,0.4 --fields 5,10,20,25,35 --anode-angle 12
```

El ángulo de 12 grados y el aluminio fijo de 3,1 mm son hipótesis de investigación de este ejemplo, NO mediciones del Azurion. Sustituirlos por la caracterización del tubo y del haz antes de interpretar resultados clínicamente. `--al` describe las familias de Al TOTAL del espectro; no sumar de nuevo el filtro de evento al generar una fila.

Los archivos ya calculados se incluyen en `skindose_verifier/data/derived_azurion_spekpy`. Permiten correr sin instalar SpekPy; solo regenerar o comprobar espectros requiere esa dependencia. La configuración de investigación es `config/equipo_philips_azurion_espectral_investigacion.json`. Utiliza el proveedor existente `tabulated_water` y selección `fixed_plus_event`: fijo 3,1 + añadido 1,0 selecciona Al total 4,1. Como el nuevo dataset incluye 0,4 mm Cu en esa familia, esa combinación deja de estar fuera de dominio. Los campos fuera de 5–35 cm, otras calidades o una geometría inválida siguen siendo limitaciones independientes. No se cambia el RDSR ni se aplican constantes históricas como fallback.

Si copias el JSON a ejemplo07, cambia sus rutas relativas, por ejemplo:
`dataset_directory: ../../DermaDose/skindose_verifier/data/derived_azurion_spekpy`
y adapta `patient.mesh_path` a tu STL. Conserva una sola configuración en la carpeta si utilizas el lanzador antiguo.

## Medidas HVL

Formato de archivo JSON: lista de objetos con `kvp`, `total_al_mm`, `cu_mm`, `hvl1_mm_al`, todos numéricos. Opcionalmente añadir procedencia/fecha para conservarlos en el archivo original. No proporcionar valores calculados haciéndolos pasar por medidas.

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier.spectral build --out tablas_contraste --kvps 80 --al 4.1 --cu 0.4 --fields 5,10,20,25,35 --anode-angle 12 --measured-hvl medidas_hvl.json --hvl-tolerance 0.05
```

Si una HVL medida difiere más de la tolerancia elegida (5 % por defecto, criterio de trabajo configurable), la generación se detiene sin publicar el dataset. Las calidades sin medida quedan marcadas, y el conjunto no figura como completamente contrastado. El prototipo compara HVL1; HVL2 calculada se registra para un contraste adicional manual. El dataset suministrado NO dispone de medidas HVL del hospital.

## Fuentes

- Benmakhlouf, Fransson y Andreo: Report KS-ASF-201101-IR, apartado 3.3, ecuaciones 11–13 y tabla 1. Informe aportado por el usuario.
- Benmakhlouf et al. PMB 2011: https://doi.org/10.1088/0031-9155/56/22/012
- Bujila et al. Validation of SpekPy (2020): https://pubmed.ncbi.nlm.nih.gov/32512239/
- Poludniowski et al. SpekPy v2.0 (2021): https://pubmed.ncbi.nlm.nih.gov/33993511/

El método mantiene el fantoma de agua, incidencia normal y corrección central del informe. No modela retrodispersión anatómica ni endurecimiento espectral de mesa/colchón; esos efectos y el registro espacial requieren evaluación independiente.

El dataset de ejemplo cubre también 46 kV, presente en un evento del RDSR aportado y ausente de la tabla clínica original; esto utiliza los nodos monoenergéticos dentro de dominio, no una extrapolación de kVp de aquella tabla clínica.
