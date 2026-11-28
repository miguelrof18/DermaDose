# Verificación con el conjunto RFM 2024

Gañán Mora A, Sánchez Casanueva RM, Fernández Soto JM. *Metodología para la validación de programas de cálculo de dosis en piel*. Rev Fis Med 2024;25(2):25–34.

El conjunto consta de 27 irradiaciones de un maniquí de agua ISO 4037 de 30 × 30 × 15 cm: los casos 1–11 en un Canon Infinix-8000 y los 12–27 en un Philips Azurion 7 M20. Para cada irradiación se dispone de su RDSR y de la dosis en agua medida con una cámara Radcal 10X6-60. La cámara se sitúa sobre la colchoneta, en el centro del campo y a la altura del IRP (IRP + 11 cm en los casos 24–27). Hay una discrepancia si la diferencia supera el 6,5 % (68 %) o el 13 % (95 %). Los RDSR son material suplementario del artículo y no se incluyen en el repositorio.

## Modelo de la medida

| Elemento | Elección | Fuente |
|---|---|---|
| Superficie | Plano de agua a la altura medida; la mesa y la colchoneta preceden a la cámara | Protocolo del artículo |
| CF | 1: el Ka,r de los RDSR ya está corregido por la calibración de la cámara de transmisión | Artículo, *Material y métodos* |
| Transmisión de mesa | Ajuste genérico de DeLorenzo et al. (2018) en función de kVp, Cu, área y ángulo; incertidumbre típica del 10 % | Parametrización de SkinDose2D, introducida como factores tabulados |
| k_BS · k_f | Canon: Benmakhlouf 2011 con 3,1 mm Al fijos y el Al del evento. Azurion: tablas `derived_azurion_hvl` (4,27 mm Al efectivos, más el Al del evento) | Perfiles de DermaDose |
| Distancias | Canon: foco–IRP 546 mm, foco–isocentro 696 mm (sus RDSR no las incluyen); área colimada referida al detector, como en el resto de equipos. Azurion: las del RDSR | Tabla 1 del artículo |
| Haces angulados | Media sobre la cara de la cámara (radio supuesto de 43,5 mm) | Criterio del artículo |

## Resultado (noviembre de 2026)

En los 27 casos la diferencia queda dentro de ±13 %. La media es de +3,0 %, con una desviación típica del 3,3 % y un máximo de −10,3 % (caso 10). Cuatro casos superan ±6,5 %: el 9, el 10 (−10,3 %), el 23 y el 27. Como referencia, los valores publicados de SkinDose2D dan −1,2 % ± 4,2 %, con cinco casos fuera de ±6,5 %.

- **Sesgo positivo, mayor con el kVp.** Al subir de 66 a 96 kV, la diferencia en el Azurion pasa del 0 al +7 %. La transmisión implícita en las medidas no crece con el kVp, mientras que el ajuste genérico de DeLorenzo crece un 3 %. Esto es compatible con la dependencia energética de la cámara (±5 %) y con la diferencia entre una mesa concreta y el ajuste genérico (10 %). La corrección consiste en medir la transmisión de la mesa y la colchoneta de la sala e introducirla como factor tabulado.
- **Casos 9 y 10 (campos del Canon con SID 910 mm).** Respecto del caso 1, la dosis medida por unidad de Ka,r varía un −11 % y un +18 %, mientras que k_BS y la transmisión de mesa predicen −7 % y +3,5 %. SkinDose2D presenta la misma discrepancia (+7 % y −9 %). Ninguna interpretación del área colimada la explica, ya que el efecto es inferior al 2 %. La causa más probable está en el Ka,r del RDSR con SID 910 mm o en las propias medidas.
- **Ángulos.** Los casos angulados 13–17 quedan entre 1,2 y 2,8 puntos por debajo del caso base del Azurion. La dependencia angular se reproduce dentro de la incertidumbre de la medida.

## Transmisión de mesa en los perfiles

La verificación anterior usa el ajuste de DeLorenzo (2018) solo en la configuración del estudio. Los perfiles del programa usan una transmisión fija, `table_normal = 0,8`, que no está medida en ninguna sala. Para saber qué supone esa hipótesis, se repiten los 27 casos con T₀ = 0,8 y la misma geometría.

| Corrección angular (`table_angular`) | Media ± DE (27 casos) | Fuera de ±6,5 % / ±13 % | Tendencia con el ángulo (casos 13–18) | Caso 14 (45°) |
|---|---|---|---|---|
| Ninguna (1, perfiles hasta 0.8.0) | −1,1 % ± 4,8 % | 4 / 1 | +2,5 % por cada 10° | +6,3 % |
| Lámina, T₀^(1/cos θ − 1) (perfiles desde 0.8.1) | −2,1 % ± 4,2 % | 3 / 1 | −0,3 % por cada 10° | −3,1 % |
| f_θ de SkinDose2D (DeLorenzo 2018) | −2,1 % ± 4,2 % | 3 / 1 | −0,3 % por cada 10° | −2,8 % |
| Fθ de DeLorenzo y Goode (2019), ec. 7 | −3,8 % ± 5,2 % | 7 / 2 | −3,8 % por cada 10° | −15,7 % |

La tendencia se ajusta sobre los casos 13–18. Con la corrección de lámina, los oblicuos 13–18 quedan entre −1,8 % y −3,4 %, y la desviación típica en el Azurion baja de 3,6 % a 2,3 %. Los 20 casos a incidencia normal no cambian. El caso 19 (+5,6 %) repite la geometría del 18 con 0,1 mm Cu en lugar de 0,4 mm. Su residuo se debe a la calidad del haz, que una T fija no recoge, y no al ángulo. El Fθ de 2019 se midió con una colchoneta de 8 cm y con parte del haz entrando por el lateral del PMMA, y sobrecorrige.

Con la corrección de lámina, la media baja a −2 % porque los oblicuos ya no compensan a los casos normales. Con T₀ = 0,82, la media del Azurion queda en +0,9 % y ningún caso supera ±13 %. Los perfiles mantienen T₀ = 0,8 hasta que se mida la mesa y la colchoneta de la sala en la geometría clínica: con el dosímetro sobre la colchoneta y retrodispersor, porque la transmisión medida en aire libre no recoge la dispersión hacia delante de la mesa (DeLorenzo y Goode, J Appl Clin Med Phys 2019;20(12):159–168).

## Cambios motivados por el estudio

- `surface_factors.aluminium_interpolation = linear_between_families`: interpola linealmente entre familias de Al adyacentes, sin extrapolar. Hace posible el caso 6 (1,8 mm Al añadidos).
- Perfil de mesa del Azurion: la columna lateral de `S` vale +1 (hacia la cabecera) según los casos 16–17. Los casos 13–15 y 18–19 confirman el signo longitudinal y los 24–27 el de la altura, con ganancia unitaria dentro de ±1 cm.
- Perfiles (`azurion_benmakhlouf_spekpy`, `ejemplo/equipo_azurion`, `benmakhlouf_tabulado`, `equipo_benmakhlouf_investigacion`): `table_angular` pasa de 1 a la tabla de lámina `SLAB-T0.80` (0–90° cada 1°, constante desde 80°). Los perfiles sintéticos de autoconsistencia conservan 1.
