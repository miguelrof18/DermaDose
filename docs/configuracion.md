# Configuración de DermaDose

## Qué archivo usar

| Configuración | Uso | Modelo dosimétrico |
|---|---|---|
| [`examples/ejemplo/equipo_azurion.json`](../examples/ejemplo/equipo_azurion.json) | Perfil de investigación para RDSR del Philips Azurion con tubo MRC 200 0407. Admite XCAT o un STL de Slicer mediante el lanzador. | `tabulated_water`. Tablas precalculadas con los datos monoenergéticos de Benmakhlouf y espectros de SpekPy. Supone 2,5 mm Al equivalentes de filtración fija y ánodo de 11°. |
| [`examples/equipo_benmakhlouf_investigacion.json`](../examples/equipo_benmakhlouf_investigacion.json) | Ejemplo y regresión con equipo simulado y superficie plana. | `benmakhlouf_2011`. Usa las tablas publicadas, para campos cuadrados y una familia de tungsteno con 3,1 mm Al asumidos. |
| [`skindose_verifier/data/equipo_sintetico.json`](../skindose_verifier/data/equipo_sintetico.json) | Perfil interno del comando `demo`. El comando puede sustituir la superficie plana por XCAT con `--model mesh`. | Factores constantes ficticios: retrodispersión 1,25 y conversión a piel 1,05. Sirve para verificar el flujo de ejecución. |

XCAT define la superficie. El JSON de equipo define la dosimetría. Con XCAT y un RDSR propio, sigue el [ejemplo Azurion](../examples/ejemplo/README.md). Para comprobar la instalación, sigue la [demo XCAT](../examples/xcat/README.md).

## Identificación y geometría

| Campo | Significado |
|---|---|
| `version` | Identificador de revisión del perfil. Permite reproducir el cálculo. |
| `manufacturer`, `model` | Fabricante y modelo de equipo que representa el perfil. |
| `validation_status` | Estado de caracterización o validación declarado. Los perfiles publicados son de investigación o simulación. |
| `dose_medium` | Medio de la dosis calculada, `water` o `skin` según el modelo. |
| `geometry.source_to_isocenter_mm`, `source_to_detector_mm`, `source_to_irp_mm` | Distancias de referencia desde el foco al isocentro, detector y punto de referencia intervencionista, respectivamente. La geometría efectiva por evento se obtiene del RDSR según las reglas del perfil. |
| `geometry.field_reference` | Plano de referencia de las dimensiones del campo, detector o isocentro. |
| `geometry.primary_sign`, `secondary_sign` | Signos de los ángulos del arco en la convención interna. |
| `geometry.rotation_order` | Orden de aplicación de los giros primario y secundario. |
| `geometry.isocenter_mm` | Coordenadas internas del isocentro. Si se omite, se toma el origen. |
| `geometry.table_coordinates` | Sistema de posiciones de mesa esperado. `legacy` usa las coordenadas convencionales del equipo e `isocenter` las relativas al isocentro. |
| `geometry.table_origin_mm`, `table_matrix` | Origen de referencia de mesa y transformación de sus ejes al sistema interno. |
| `geometry.table_plane_z_mm` | Altura del plano de la mesa en el modelo de referencia. |
| `geometry.table_half_width_mm`, `table_half_length_mm` | Semianchura y semilongitud del tablero para detectar si un rayo lo atraviesa. |

## Superficie y factores de dosis

| Campo | Significado |
|---|---|
| `patient.type` | Tipo de superficie. Los perfiles sintéticos parten de un plano y el lanzador permite seleccionar una malla. |
| `patient.width_mm`, `length_mm`, `z_mm`, `normal` | Dimensiones, altura y normal del plano sintético. |
| `patient.thickness_mm`, `resolution_mm` | Espesor del modelo y resolución nominal del muestreo. |
| `patient.origin_on_table_mm` | Desplazamiento de la superficie respecto a la mesa. La malla puede tener además su propia transformación. |
| `factors.cf` | Factor de corrección del kerma registrado. |
| `factors.table_normal`, `table_angular` | Transmisión normal de mesa y corrección angular configurada, aplicables a rayos que atraviesan el tablero. |
| `factors.bsf`, `f_skin` | Retrodispersión y conversión a piel explícitas en la demo sintética. |
| `value`, `provenance`, `identifier`, `relative_standard_uncertainty` | Valor del factor, procedencia, identificador e incertidumbre estándar relativa. `null` indica que no se ha declarado la incertidumbre. |
| `tube_output.reference_distance_mm`, `yield_gy_per_mas` | Distancia y rendimiento de referencia, en Gy/mAs, para la estimación alternativa de kerma por rendimiento. |
| `rdsr_relative_standard_uncertainty` | Incertidumbre estándar relativa asignada al dato dosimétrico del RDSR. |
| `allow_uniform_sweep` | Permite aproximar un barrido angular con reparto uniforme cuando el cálculo lo admite. |

## Tablas y perfil Azurion

| Campo | Significado |
|---|---|
| `surface_factors.provider` | Selecciona tablas publicadas o derivadas. El [metadato del Azurion](../skindose_verifier/data/derived_mrc200_nominal/metadata.json) documenta su derivación con SpekPy y la ecuación 13 de Benmakhlouf. |
| `surface_factors.dataset_directory` | Carpeta de la tabla derivada y su metadato. La ruta relativa se interpreta desde el JSON de equipo. |
| `surface_factors.inherent_al_mm`, `target` | Filtración fija efectiva en mm Al equivalentes y material del blanco. En el perfil Azurion, el valor representa el conjunto nominal. |
| `surface_factors.aluminium_policy` | Regla para combinar el aluminio fijo con el aluminio registrado por evento. |
| `surface_factors.expected_anode_angle_deg` | Ángulo de ánodo exigido para la tabla seleccionada. |
| `surface_factors.field_shape_policy`, `out_of_domain` | Regla para representar la forma del campo y comportamiento si los parámetros quedan fuera del dominio tabulado. `equal_area` usa un cuadrado de área equivalente. |
| `surface_factors.incident_kerma_basis`, `quality_basis` | Definición del kerma incidente y aproximación de calidad de haz usada para consultar la tabla. |
| `surface_factors.provenance`, `allow_water_phantom_approximation` | Origen e hipótesis del modelo y permiso explícito para usar la aproximación de fantoma de agua. |
| `azurion_profile.name`, `tube_model` | Identifican las reglas Philips y el conjunto de tubo asociado. |
| `azurion_profile.allow_research_assumptions` | Permite las hipótesis del perfil de investigación. |
| `azurion_profile.fixed_filtration_characterized`, `geometry_characterized`, `calibration_characterized` | Declaran si se han caracterizado la filtración fija, la geometría y la calibración del equipo. |
| `azurion_profile.reference_table_from_first_event` | Usa la posición de mesa del primer evento como origen provisional del caso. |
| `tube_characterization.model`, `fixed_effective_al_mm`, `filtration_basis` | Identidad del tubo, aluminio fijo efectivo y base de ese valor. Deben ser coherentes con la tabla espectral seleccionada. |

Mantén los parámetros de sala en el perfil de equipo y la transformación específica del paciente en los datos del caso. Antes de interpretar la colocación como representativa de una intervención, comprueba con un fantoma la lateralidad, los signos angulares y los desplazamientos de mesa.
