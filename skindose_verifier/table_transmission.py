"""Generic table + mattress transmission as tabulated DermaDose factors.

DeLorenzo MC, Yang K, Li X, Liu B. Comprehensive evaluation of broad-beam
transmission of patient supports. Med Phys 2018;45(4). doi:10.1002/mp.12803.
The parametrisation is the one coded in SkinDose2D (TableTransmision.m, 2023),
used by Gañán Mora et al. (Rev Fis Med 2024;25(2):25-34):

    T = f_kvf(kVp, Cu) * f_tc(area) * f_ang(theta)

It is a generic fit over several tables, not a measurement of a given table:
the default relative standard uncertainty (0.10) is the one of that article
(Table 4). Replace it with measured factors of the installation when available.

    python -m skindose_verifier.table_transmission --salida factores_mesa.json
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

SOURCE = ("DeLorenzo et al., Med Phys 2018;45(4), doi:10.1002/mp.12803; parametrisation of "
          "SkinDose2D TableTransmision.m (Gañán Mora et al., Rev Fis Med 2024). Generic fit, not a "
          "measurement of this table.")
AREA_CAP_CM2 = 580.  # f_tc is held constant above 0.058 m2
ANGLE_LIMIT_DEG = 80.  # f_ang is held at its 80 degree value beyond (SkinDose2D jumps to 1)

KVP_AXIS = [float(k) for k in range(40, 152, 2)]
CU_AXIS = [0., .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .6, .7, .8, .9, 1.]
AREA_AXIS = [1.] + [float(a) for a in range(10, 590, 10)] + [600., 10000.]
ANGLE_AXIS = [float(a) for a in range(0, 81)] + [90.]


def kvp_filtration_term(kvp: float, cu_mm: float) -> float:
    if not (math.isfinite(kvp) and kvp > 0 and math.isfinite(cu_mm) and cu_mm >= 0):
        raise ValueError("kVp y Cu deben ser finitos, kVp > 0 y Cu >= 0.")
    return 1 + 0.1321*math.log(1.5892*kvp + 0.1321) + 5.043e-2*math.log(1.5892*cu_mm + 5.043e-2)


def field_area_term(area_cm2: float) -> float:
    if not (math.isfinite(area_cm2) and area_cm2 > 0):
        raise ValueError("Área de campo no positiva.")
    tc = min(area_cm2, AREA_CAP_CM2)/1e4
    return -12.47*tc**2 + 1.7081*tc + 0.4890


def angular_term(theta_deg: float) -> float:
    if not (math.isfinite(theta_deg) and 0 <= theta_deg <= 90):
        raise ValueError("Ángulo de incidencia en la mesa fuera de [0, 90] grados.")
    a = min(theta_deg, ANGLE_LIMIT_DEG)
    return -4.90e-05*a**2 + 3.07e-04*a + 1.0


def transmission(kvp: float, cu_mm: float, area_cm2: float, theta_deg: float = 0.) -> float:
    return kvp_filtration_term(kvp, cu_mm)*field_area_term(area_cm2)*angular_term(theta_deg)


def delorenzo_factor_blocks(relative_standard_uncertainty: float = 0.10) -> dict:
    """`table_normal` (kVp, Cu, area at the table) and `table_angular` (incidence) blocks."""
    normal = [[[kvp_filtration_term(k, c)*field_area_term(a) for a in AREA_AXIS] for c in CU_AXIS] for k in KVP_AXIS]
    return {
        "table_normal": {
            "axes": {"kvp": KVP_AXIS, "filtration_mm_cu": CU_AXIS, "field_area_cm2": AREA_AXIS},
            "grid": normal, "provenance": SOURCE, "identifier": "DELORENZO2018-SKINDOSE2D-NORMAL",
            "relative_standard_uncertainty": relative_standard_uncertainty},
        "table_angular": {
            "axes": {"table_angle_deg": ANGLE_AXIS}, "grid": [angular_term(a) for a in ANGLE_AXIS],
            "provenance": SOURCE + " Angular term held at its 80 degree value beyond 80 degrees.",
            "identifier": "DELORENZO2018-SKINDOSE2D-ANGULAR", "relative_standard_uncertainty": None},
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Factores genéricos de transmisión de mesa (DeLorenzo 2018).")
    parser.add_argument("--salida", type=Path, help="JSON con los bloques table_normal y table_angular.")
    parser.add_argument("--incertidumbre", type=float, default=0.10, help="Incertidumbre relativa típica (0.10).")
    args = parser.parse_args(argv)
    text = json.dumps(delorenzo_factor_blocks(args.incertidumbre), indent=1)
    if args.salida:
        args.salida.write_text(text, encoding="utf-8")
        print(f"Bloques de factores escritos en {args.salida}. Cópielos en 'factors' de la configuración.")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
