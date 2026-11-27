# Verificación con el conjunto RFM 2024

Gañán Mora A, Sánchez Casanueva RM, Fernández Soto JM. *Metodología para la validación de programas de cálculo de dosis en piel*. Rev Fis Med 2024;25(2):25–34.

El conjunto consta de 27 irradiaciones de un maniquí de agua ISO 4037 de 30 × 30 × 15 cm: los casos 1–11 en un Canon Infinix-8000 y los 12–27 en un Philips Azurion 7 M20. Para cada irradiación se dispone de su RDSR y de la dosis en agua medida con una cámara Radcal 10X6-60. La cámara se sitúa sobre la colchoneta, en el centro del campo y a la altura del IRP (IRP + 11 cm en los casos 24–27). Hay una discrepancia si la diferencia supera el 6,5 % (68 %) o el 13 % (95 %). Los RDSR son material suplementario del artículo y no se incluyen en el repositorio.

```powershell
python -m skindose_verifier.verification_rfm2024 SR_Finals_240501 --salida resultados_rfm2024
python -m skindose_verifier.verification_rfm2024 SR_Finals_240501 --altura rdsr   # altura de las lecturas de mesa
```

El módulo escribe `resultados_rfm2024.csv` y `.json` con la dosis de DermaDose, la medida, el valor publicado de SkinDose2D y el veredicto de cada caso.

## Modelo de la medida

| Elemento | Elección | Fuente |
|---|---|---|
| Superficie | Plano de agua a la altura medida; la mesa y la colchoneta preceden a la cámara | Protocolo del artículo |
| CF | 1: el Ka,r de los RDSR ya está corregido por la calibración de la cámara de transmisión | Artículo, *Material y métodos* |
| Transmisión de mesa | Ajuste genérico de DeLorenzo et al. (2018) en función de kVp, Cu, área y ángulo; incertidumbre típica del 10 % | `table_transmission.py` |
| k_BS · k_f | Canon: Benmakhlouf 2011 con 3,1 mm Al fijos y el Al del evento. Azurion: tablas `derived_azurion_hvl` (4,27 mm Al efectivos, más el Al del evento) | Perfiles de DermaDose |
| Distancias | Canon: foco–IRP 546 mm, foco–isocentro 696 mm (sus RDSR no las incluyen); área colimada referida a 1000 mm. Azurion: las del RDSR | Tabla 1 del artículo y casos 9–10 |
| Haces angulados | Media sobre la cara de la cámara (radio supuesto de 43,5 mm) | Criterio del artículo |

## Resultado (noviembre de 2026)

En los 27 casos la diferencia queda dentro de ±13 %. La media es de +2,9 %, con una desviación típica del 3,4 % y un máximo de −10,3 % (caso 10). Tres casos superan ±6,5 %: el 10 (−10,3 %), el 23 (+6,8 %) y el 27 (+6,7 %). Con la altura de las lecturas de mesa (`--altura rdsr`), la media baja a +2,4 % y solo el caso 10 supera ±6,5 %. Como referencia, los valores publicados de SkinDose2D dan −1,2 % ± 4,2 %, con cinco casos fuera de ±6,5 %.

- **Sesgo positivo, mayor con el kVp.** Al subir de 66 a 96 kV, la diferencia en el Azurion pasa del 0 al +7 %. La transmisión implícita en las medidas no crece con el kVp, mientras que el ajuste genérico de DeLorenzo crece un 3 %. Esto es compatible con la dependencia energética de la cámara (±5 %) y con la diferencia entre una mesa concreta y el ajuste genérico (10 %). La corrección consiste en medir la transmisión de la mesa y la colchoneta de la sala e introducirla como factor tabulado.
- **Casos 9 y 10 (campos del Canon con SID 910 mm).** Respecto del caso 1, la dosis medida por unidad de Ka,r varía un −11 % y un +18 %, mientras que k_BS y la transmisión de mesa predicen −7 % y +3,5 %. SkinDose2D presenta la misma discrepancia (+7 % y −9 %). Ninguna interpretación del área colimada la explica, ya que el efecto es inferior al 2 %. La causa más probable está en el Ka,r del RDSR con SID 910 mm o en las propias medidas.
- **Ángulos.** Los casos angulados 13–17 quedan entre 1,2 y 2,8 puntos por debajo del caso base del Azurion. La dependencia angular se reproduce dentro de la incertidumbre de la medida.

## Cambios motivados por el estudio

- `surface_factors.aluminium_interpolation = linear_between_families`: interpola linealmente entre familias de Al adyacentes, sin extrapolar. Hace posible el caso 6 (1,8 mm Al añadidos).
- `geometry.collimated_area_reference_mm`: indica la distancia a la que el equipo refiere el área colimada cuando no proporciona anchura ni altura (Canon).
- `table_transmission.py`: genera bloques de factores con el ajuste de DeLorenzo y los términos documentados.
- Perfil de mesa del Azurion: la columna lateral de `S` vale +1 (hacia la cabecera) según los casos 16–17. Los casos 13–15 y 18–19 confirman el signo longitudinal y los 24–27 el de la altura, con ganancia unitaria dentro de ±1 cm.
