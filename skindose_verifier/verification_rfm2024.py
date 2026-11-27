"""Verification against the RFM 2024 reference set (Gañán Mora, Sánchez Casanueva, Fernández Soto).

"Metodología para la validación de programas de cálculo de dosis en piel",
Rev Fis Med 2024;25(2):25-34. 27 irradiations of a 30x30x15 cm water phantom
(Canon Infinix-8000, cases 1-11; Philips Azurion 7 M20, cases 12-27) with the
RDSR of each irradiation and the dose to water measured with a Radcal 10X6-60
at the centre of the field, at the IRP height (IRP + 11 cm in cases 24-27).

The RDSR files are the article's supplementary material and are not shipped.

    python -m skindose_verifier.verification_rfm2024 SR_Finals_240501 --salida resultados_rfm2024

Model of the measurement: plane water surface at the measured height, beam
passing through table and mattress, CF = 1 (the RDSR Ka,r are already corrected
by the transmission-chamber calibration), DeLorenzo table transmission, the
equipment's surface factors, and for angled beams the mean over the chamber face.
Tolerances: |difference| > 6.5 % (68 %) or > 13 % (95 %) is a discrepancy.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import replace
from importlib.resources import files
from pathlib import Path

import numpy as np

from .config import Equipment
from .domain import InputError
from .engine import calculate, normalize
from .geometry import beam_from_event
from .rdsr import read_rdsr
from .table_transmission import delorenzo_factor_blocks

# Case numbering of the RDSR file names (CASE_n). Table 5 of the article lists
# cases 8-10 in another order; the measured dose identifies each one unambiguously.
# case: (equipment, purpose, Dw mGy, 95 % interval, kf, SkinDose2D mGy, plane above IRP mm)
REFERENCE = {
    1: ("Canon", "Base Canon (80 kV, 0,2 Cu, 16,4 cm)", 13.26, (12.20, 14.39), 1.032, 13.3, 0.),
    2: ("Canon", "Filtración 0,3 Cu", 13.39, (12.32, 14.53), 1.034, 13.4, 0.),
    3: ("Canon", "Filtración 0,5 Cu", 13.78, (12.67, 14.95), 1.037, 13.7, 0.),
    4: ("Canon", "60 kV", 12.09, (11.12, 13.11), 1.022, 11.6, 0.),
    5: ("Canon", "100 kV", 14.01, (12.89, 15.20), 1.041, 14.5, 0.),
    6: ("Canon", "60 kV, 1,8 mm Al", 11.69, (10.75, 12.68), 1.020, 11.8, 0.),
    7: ("Canon", "120 kV, 0,5 Cu", 17.17, (15.80, 18.63), 1.055, 18.5, 0.),
    8: ("Canon", "Espesor +4 cm PMMA", 13.02, (11.98, 14.13), 1.032, 13.1, 0.),
    9: ("Canon", "Campo 11×11 cm", 12.44, (11.44, 13.49), 1.033, 13.3, 0.),
    10: ("Canon", "Campo 16,4×22 cm", 15.17, (13.96, 16.46), 1.032, 13.9, 0.),
    11: ("Canon", "Espesor −6 cm de agua", 14.95, (13.75, 16.22), 1.032, 14.9, 0.),
    12: ("Philips", "Base Philips (FP2, 21,7×16,8 cm)", 11.52, (10.63, 12.45), 1.026, 12.0, 0.),
    13: ("Philips", "Primario 30°", 10.15, (9.35, 10.99), 1.028, 9.6, 0.),
    14: ("Philips", "Primario 45°", 10.33, (9.51, 11.18), 1.029, 10.4, 0.),
    15: ("Philips", "Primario −15°", 9.90, (9.12, 10.72), 1.027, 9.4, 0.),
    16: ("Philips", "(−15°, 15°)", 10.00, (9.21, 10.83), 1.027, 9.5, 0.),
    17: ("Philips", "(−15°, 30°)", 9.59, (8.84, 10.39), 1.028, 9.1, 0.),
    18: ("Philips", "−30°, campo 12,8 cm, FP2", 10.22, (9.42, 11.07), 1.030, 10.0, 0.),
    19: ("Philips", "−30°, campo 12,8 cm, FP1", 10.36, (9.55, 11.22), 1.029, 10.1, 0.),
    20: ("Philips", "66 kV, campo 14,8 cm", 10.44, (9.64, 11.29), 1.027, 9.6, 0.),
    21: ("Philips", "76 kV", 10.45, (9.65, 11.30), 1.033, 10.1, 0.),
    22: ("Philips", "86 kV", 14.29, (13.19, 15.45), 1.038, 14.4, 0.),
    23: ("Philips", "96 kV", 14.07, (12.99, 15.21), 1.044, 14.2, 0.),
    24: ("Philips", "IRP + 11 cm, 66 kV", 7.92, (7.33, 8.55), 1.027, 7.3, 110.),
    25: ("Philips", "IRP + 11 cm, 76 kV", 8.80, (8.14, 9.49), 1.033, 8.5, 110.),
    26: ("Philips", "IRP + 11 cm, 85 kV", 10.66, (9.86, 11.50), 1.038, 10.7, 110.),
    27: ("Philips", "IRP + 11 cm, 95 kV", 12.25, (11.33, 13.22), 1.043, 12.4, 110.),
}
TOLERANCE_68, TOLERANCE_95 = 6.5, 13.0
BASE_CASE = {"Canon": 1, "Philips": 12}
# Effective radius of the Radcal 10X6-60 face: not given in the article; assumption.
CHAMBER_RADIUS_MM = 43.5

# Table 1 of the article: focus-IRP and focus-detector distances; IRP 15 cm before the isocentre.
GEOMETRY = {
    "Canon": {"source_to_isocenter_mm": 696., "source_to_detector_mm": 1000., "source_to_irp_mm": 546.,
              # The Canon RDSR reports the collimated area at 1 m from the focus, not at the
              # detector: cases 9-10 (SID 910 mm) reproduce Table 2 only with this reference.
              "collimated_area_reference_mm": 1000.},
    "Philips": {"source_to_isocenter_mm": 810., "source_to_detector_mm": 1181.4, "source_to_irp_mm": 660.},
}


def _surface(device: str) -> dict:
    common = {"target": "W", "field_shape_policy": "equal_area", "out_of_domain": "error",
              "incident_kerma_basis": "free_in_air_without_backscatter",
              "quality_basis": "event_kvp_filtration_proxy_at_surface",
              "allow_water_phantom_approximation": True, "aluminium_policy": "fixed_plus_event",
              "aluminium_interpolation": "linear_between_families"}
    if device == "Canon":
        return {**common, "provider": "benmakhlouf_2011", "inherent_al_mm": 3.1,
                "provenance": "Benmakhlouf 2011 Table 3; Canon fixed filtration not characterised (3.1 mm Al assumed)."}
    return {**common, "provider": "tabulated_water", "inherent_al_mm": 4.270081,
            "dataset_directory": str(files("skindose_verifier").joinpath("data/derived_azurion_hvl")),
            "provenance": "DermaDose Azurion HVL dataset (4.27 mm Al effective) + 1 mm Al of the event."}


def study_configuration(device: str, plane_z_mm: float, *, surface: dict | None = None,
                        table_factors: dict | None = None) -> dict:
    """Equipment configuration that models the RFM 2024 set-up for one equipment."""
    if device not in GEOMETRY:
        raise InputError(f"Equipo desconocido en el conjunto RFM 2024: {device}.")
    geometry = {**GEOMETRY[device], "field_reference": "detector", "primary_sign": 1, "secondary_sign": 1,
                "rotation_order": "primary_then_secondary", "isocenter_mm": [0, 0, 0],
                "table_coordinates": "isocenter", "table_origin_mm": [0, 0, 0],
                "table_matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                # Chamber on the mattress: the whole field crosses table and mattress first.
                "table_plane_z_mm": plane_z_mm - 0.5, "table_half_width_mm": 400, "table_half_length_mm": 1000}
    return {
        "version": f"RFM2024-{device}", "manufacturer": device, "model": "RFM2024 verification set",
        "validation_status": "verification_set_Rev_Fis_Med_2024", "dose_medium": "water", "geometry": geometry,
        "patient": {"type": "plane", "width_mm": 700, "length_mm": 700, "z_mm": plane_z_mm, "normal": [0, 0, -1],
                    "thickness_mm": 150, "origin_on_table_mm": [0, 0, 0], "resolution_mm": 2.5},
        "factors": {"cf": {"value": 1.0, "identifier": "RFM2024-CF1", "relative_standard_uncertainty": 0.037,
                           "provenance": "RFM 2024: RDSR Ka,r already corrected by the transmission-chamber calibration."},
                    **(table_factors or delorenzo_factor_blocks())},
        "rdsr_relative_standard_uncertainty": None, "allow_uniform_sweep": False,
        "surface_factors": surface or _surface(device),
    }


def case_from_name(path: Path) -> int:
    match = re.search(r"CASE_(\d+)", path.name)
    if not match or int(match.group(1)) not in REFERENCE:
        raise InputError(f"{path.name}: no corresponde a un caso del conjunto RFM 2024.")
    case = int(match.group(1))
    measured = re.search(r"MeasDose=([\d.]+)", path.name)
    if measured and not math.isclose(float(measured.group(1)), REFERENCE[case][2], abs_tol=0.006):
        raise InputError(f"{path.name}: la dosis del nombre no coincide con la tabla de referencia.")
    return case


def plane_height(case: int, events, heights: str, base_height_mm: float | None) -> float:
    """Height of the phantom entrance surface (engine z, isocentre at 0)."""
    device, *_, above_irp = REFERENCE[case]
    irp = -(GEOMETRY[device]["source_to_isocenter_mm"] - GEOMETRY[device]["source_to_irp_mm"])
    if heights == "protocolo" or base_height_mm is None:
        return irp + above_irp
    if heights != "rdsr":
        raise InputError("--altura debe ser protocolo o rdsr.")
    readings = {float(e.table_mm[2]) for e in events if e.table_system == "legacy" and e.table_mm}
    if len(readings) != 1:
        raise InputError(f"Caso {case}: altura de mesa ausente o variable entre eventos.")
    # Height reading positive upwards (verified with cases 24-27), relative to the base case at the IRP.
    return irp + readings.pop() - base_height_mm


def evaluate_case(path: Path, *, heights: str = "protocolo", base_height_mm: float | None = None,
                  chamber_radius_mm: float = CHAMBER_RADIUS_MM, configuration=None) -> dict:
    case = case_from_name(path)
    device, purpose, dw, interval, kf, skindose2d, _ = REFERENCE[case]
    data = read_rdsr(path)
    z = plane_height(case, data.events, heights, base_height_mm)
    cfg = (configuration or study_configuration)(device, z)
    equipment = Equipment.from_dict(cfg, base_directory=Path("."))
    # The chamber moves with the table to the field centre: a horizontal shift of the
    # plane phantom does not change the dose there, so the table is placed at the origin.
    events = [replace(e, table_mm=(0., 0., 0.), table_system="isocenter") for e in data.events]
    result = calculate(events, equipment, "skin_dose", detailed=False)
    points = result.surface.points_mm
    dose, details = 0., []
    for index, event in enumerate(events):
        beam = beam_from_event(normalize(event, equipment), equipment)
        centre = beam.source_mm + (z - beam.source_mm[2])/beam.direction[2]*beam.direction
        angled = abs(event.primary_deg) > .5 or abs(event.secondary_deg) > .5
        radius = max(chamber_radius_mm if angled else 0., 3.)
        selected = np.hypot(points[:, 0]-centre[0], points[:, 1]-centre[1]) <= radius
        if not np.any(selected):
            raise InputError(f"Caso {case}: el centro del campo queda fuera del maniquí modelado.")
        value = float(np.mean(result.event_maps_gy[index][selected]))*1000.
        dose += value
        details.append({"uid": event.uid, "dose_mgy": value, "ka_r_mgy": event.ka_r_gy*1000., "kvp": event.kvp,
                        "cu_mm": event.filtration_mm_cu, "primary_deg": event.primary_deg,
                        "secondary_deg": event.secondary_deg, "field_centre_mm": centre.tolist(),
                        "chamber_radius_mm": radius})
    difference = 100.*(dose/dw - 1.)
    return {"case": case, "equipment": device, "purpose": purpose, "plane_z_mm": z, "dermadose_mgy": dose,
            "measured_dw_mgy": dw, "interval_95_mgy": list(interval), "kf_article": kf,
            "skindose2d_mgy": skindose2d, "difference_percent": difference,
            "skindose2d_difference_percent": 100.*(skindose2d/dw - 1.), "verdict": verdict(difference),
            "events": details, "rdsr_sha256": data.sha256}


def verdict(difference_percent: float) -> str:
    d = abs(difference_percent)
    if d > TOLERANCE_95:
        return "discrepancia significativa (95 %)"
    if d > TOLERANCE_68:
        return "discrepancia al 68 %, no significativa al 95 %"
    return "sin discrepancia"


def summary(rows: list[dict]) -> dict:
    d = np.array([r["difference_percent"] for r in rows])
    out = {"cases": len(rows), "mean_percent": float(d.mean()), "sd_percent": float(d.std(ddof=1)) if len(d) > 1 else 0.,
           "rms_percent": float(np.sqrt(np.mean(d**2))), "max_abs_percent": float(np.max(np.abs(d))),
           "beyond_68": int(np.sum(np.abs(d) > TOLERANCE_68)), "beyond_95": int(np.sum(np.abs(d) > TOLERANCE_95))}
    for device in ("Canon", "Philips"):
        sub = [r["difference_percent"] for r in rows if r["equipment"] == device]
        if sub:
            out[f"mean_percent_{device}"] = float(np.mean(sub))
    return out


def run(folder: Path, *, heights: str = "protocolo", chamber_radius_mm: float = CHAMBER_RADIUS_MM) -> list[dict]:
    paths = sorted((p for p in Path(folder).rglob("*") if p.is_file() and "CASE_" in p.name), key=case_from_name)
    if not paths:
        raise InputError(f"No se encontraron RDSR CASE_n en {folder}.")
    base = {}
    if heights == "rdsr":
        for p in paths:
            case = case_from_name(p)
            if case in BASE_CASE.values():
                base[REFERENCE[case][0]] = float(read_rdsr(p).events[0].table_mm[2])
    return [evaluate_case(p, heights=heights, base_height_mm=base.get(REFERENCE[case_from_name(p)][0]),
                          chamber_radius_mm=chamber_radius_mm) for p in paths]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Verificación con el conjunto RFM 2024 (27 irradiaciones).")
    parser.add_argument("carpeta", type=Path, help="Carpeta SR_Finals con los RDSR CASE_n.")
    parser.add_argument("--salida", type=Path, default=Path("resultados_rfm2024"), help="Prefijo de CSV y JSON.")
    parser.add_argument("--altura", choices=("protocolo", "rdsr"), default="protocolo",
                        help="Altura del plano: la del protocolo, o la lectura de mesa relativa al caso base.")
    parser.add_argument("--radio-camara-mm", type=float, default=CHAMBER_RADIUS_MM)
    args = parser.parse_args(argv)
    rows = run(args.carpeta, heights=args.altura, chamber_radius_mm=args.radio_camara_mm)
    stats = summary(rows)
    args.salida.with_suffix(".json").write_text(json.dumps({"summary": stats, "cases": rows}, indent=1), encoding="utf-8")
    with args.salida.with_suffix(".csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["caso", "equipo", "verificacion", "Dw_medida_mGy", "DermaDose_mGy", "diferencia_%",
                         "SkinDose2D_mGy", "SkinDose2D_%", "veredicto"])
        for r in rows:
            writer.writerow([r["case"], r["equipment"], r["purpose"], f"{r['measured_dw_mgy']:.2f}",
                             f"{r['dermadose_mgy']:.2f}", f"{r['difference_percent']:+.1f}", f"{r['skindose2d_mgy']:.1f}",
                             f"{r['skindose2d_difference_percent']:+.1f}", r["verdict"]])
    for r in rows:
        print(f"{r['case']:2d} {r['equipment']:7s} {r['purpose'][:34]:34s} Dw {r['measured_dw_mgy']:6.2f}  "
              f"DermaDose {r['dermadose_mgy']:6.2f}  {r['difference_percent']:+5.1f} %  {r['verdict']}")
    print(f"Media {stats['mean_percent']:+.1f} %, DE {stats['sd_percent']:.1f} %, máx. |dif.| {stats['max_abs_percent']:.1f} %, "
          f"fuera de ±6,5 %: {stats['beyond_68']}, fuera de ±13 %: {stats['beyond_95']}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
