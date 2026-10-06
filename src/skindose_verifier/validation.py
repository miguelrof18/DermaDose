"""Comparación de referencias experimentales (primaria) o comerciales (secundaria)."""
from __future__ import annotations

from pathlib import Path
import csv
import json
import math

import numpy as np
from scipy.spatial import cKDTree

from .config import Equipment
from .domain import InputError
from .engine import calculate
from .rdsr import read_rdsr


def compare_with_radimetrics(ours: dict, reference_kar_gy: float,
                             reference_psd_gy: float, reference_medium: str) -> dict:
    if ours["stage"] != "skin_dose":
        raise InputError("PSD comercial solo es comparable con cálculo de dosis absorbida.")
    if ours["quantity"] != f"entrance_absorbed_dose_in_{reference_medium}_Gy":
        raise InputError("La magnitud/medio dosimétrico difiere de Radimetrics.")
    if reference_kar_gy <= 0 or reference_psd_gy <= 0:
        raise InputError("Valores de referencia no positivos.")
    x, y = ours["total_rdsr_ka_r_gy"], ours["psd_point_gy"]
    return {"reference_is_not_ground_truth": True,
            "delta_kar_gy": x-reference_kar_gy, "relative_kar": x/reference_kar_gy-1,
            "delta_psd_gy": y-reference_psd_gy, "relative_psd": y/reference_psd_gy-1,
            "note": "Ka,r de este campo procede del mismo RDSR; usar la rama tube_output para independencia."}


def validate_measurements(csv_path: str | Path, stage: str) -> dict:
    """CSV: test_id,rdsr_path,config_path,measured_max_gy,standard_uncertainty_gy.

    Measured quantity and measurement footprint must match requested stage and
    point-PSD definition; this function does not make that metrological decision.
    """
    path = Path(csv_path)
    with path.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise InputError("CSV experimental sin filas.")
    results = []
    for row in rows:
        if not all(row.get(k) for k in ("test_id", "rdsr_path", "config_path", "measured_max_gy", "standard_uncertainty_gy")):
            raise InputError("CSV experimental requiere test_id,rdsr_path,config_path,measured_max_gy,standard_uncertainty_gy.")
        equipment = Equipment.load(path.parent/row["config_path"])
        rdsr = read_rdsr(path.parent/row["rdsr_path"])
        equipment.require_matching_device(rdsr.manufacturer, rdsr.model)
        prediction = calculate(rdsr.events, equipment, stage).summary["psd_point_gy"]
        measured, u = float(row["measured_max_gy"]), float(row["standard_uncertainty_gy"])
        if not all(math.isfinite(v) for v in (measured, u)) or measured <= 0 or u <= 0:
            raise InputError("Medida experimental o incertidumbre inválida.")
        results.append({"test_id": row["test_id"], "predicted_gy": prediction,
                        "measured_gy": measured, "measured_standard_uncertainty_gy": u,
                        "relative_difference": prediction/measured-1,
                        "difference_over_measurement_u": (prediction-measured)/u,
                        "measurement_footprint_not_verified_by_software": True})
    differences = np.asarray([r["relative_difference"] for r in results])
    return {"cases": results, "mean_relative_difference": float(np.mean(differences)),
            "sample_std_relative_difference": None if len(results)<2 else float(np.std(differences, ddof=1)),
            "validated": False,
            "warning": "Comprobar medio, geometría, área del dosímetro, trazabilidad e incertidumbre del modelo antes de aceptar."}


def validate_map(results_path: str | Path, measurements_path: str | Path,
                 threshold_fraction: float = .5, max_distance_mm: float = 5.) -> dict:
    """Compare registered phantom point samples; alignment is supplied, never estimated here.

    CSV: x_mm,y_mm,z_mm,measured_gy,standard_uncertainty_gy,quantity.
    Both arrays are evaluated only at the supplied measurement locations.
    """
    if not (0 < threshold_fraction < 1) or not (max_distance_mm > 0 and math.isfinite(max_distance_mm)):
        raise InputError("Umbral o tolerancia espacial inválidos.")
    report_path = Path(results_path)
    result = json.loads(report_path.read_text(encoding="utf-8"))
    with np.load(report_path.parent/"mapa_dosis.npz", allow_pickle=False) as contents:
        model_points = np.asarray(contents["patient_points_mm"], dtype=float)
        doses = np.asarray(contents["dose_gy"], dtype=float).ravel()
    if len(doses) != len(model_points) or len(doses) == 0:
        raise InputError("Malla de resultados inconsistente.")
    with Path(measurements_path).open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    fields = ("x_mm", "y_mm", "z_mm", "measured_gy", "standard_uncertainty_gy", "quantity")
    if len(rows) < 4 or not all(all(row.get(key) for key in fields) for row in rows):
        raise InputError("Mapa experimental: al menos 4 puntos y columnas x_mm,y_mm,z_mm,measured_gy,standard_uncertainty_gy,quantity.")
    if any(row["quantity"] != result["quantity"] for row in rows):
        raise InputError("Magnitud o medio del mapa experimental no coincide con el resultado.")
    try:
        points = np.array([[float(row[k]) for k in ("x_mm", "y_mm", "z_mm")] for row in rows])
        measured = np.array([float(row["measured_gy"]) for row in rows])
        u = np.array([float(row["standard_uncertainty_gy"]) for row in rows])
    except ValueError as exc:
        raise InputError("Mapa experimental contiene números inválidos.") from exc
    if not (np.all(np.isfinite(points)) and np.all(np.isfinite(measured)) and np.all(np.isfinite(u))):
        raise InputError("Mapa experimental contiene datos no finitos.")
    if np.any(measured < 0) or np.any(u <= 0) or not np.any(measured > 0):
        raise InputError("Dosis experimentales no positivas o incertidumbres inválidas.")
    distances, indices = cKDTree(model_points).query(points)
    if np.any(distances > max_distance_mm) or len(np.unique(indices)) != len(indices):
        raise InputError("Registro espacial insuficiente: distancia excesiva o puntos duplicados de la malla.")
    predicted = doses[indices]
    diff = predicted-measured
    threshold = threshold_fraction * float(np.max(measured))
    high_pred = predicted >= threshold
    high_measured = measured >= threshold
    intersection = int(np.count_nonzero(high_pred & high_measured))
    denom = int(np.count_nonzero(high_pred) + np.count_nonzero(high_measured))
    idx_pred = int(np.argmax(predicted))
    idx_measured = int(np.argmax(measured))
    return {
        "quantity": result["quantity"], "sample_count": len(rows),
        "maximum_registration_distance_mm": float(np.max(distances)),
        "mean_difference_gy": float(np.mean(diff)),
        "rmse_gy": float(np.sqrt(np.mean(diff**2))),
        "peak_measured_gy": float(measured[idx_measured]),
        "peak_predicted_at_measured_locations_gy": float(predicted[idx_pred]),
        "peak_location_separation_mm": float(np.linalg.norm(points[idx_pred]-points[idx_measured])),
        "fractional_isodose_threshold_gy": threshold,
        "isodose_dice_at_sampled_points": (2*intersection/denom if denom else None),
        "mean_difference_over_measurement_u": float(np.mean(diff/u)),
        "validated": False,
        "limitations": ["Registro de puntos y magnitud dosimétrica proporcionados por el usuario",
                        "Solo puntos muestreados, no área completa de isodosis",
                        "Huella, calibración e incertidumbre del modelo no evaluadas"],
    }
