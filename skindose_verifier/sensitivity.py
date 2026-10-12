"""Explicit scenario analysis for uncertain patient/table alignment."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
import math
from typing import Sequence

import numpy as np

from .config import Equipment
from .domain import InputError, IrradiationEvent
from .engine import calculate
from .patient import sample_patient


def position_sensitivity(events: Sequence[IrradiationEvent], equipment: Equipment,
                         stage: str, offset_mm: float = 5.) -> dict:
    """Recalculate at six explicit offsets; the envelope is not a statistical CI."""
    if not math.isfinite(offset_mm) or offset_mm <= 0 or offset_mm > 100:
        raise InputError("El desplazamiento de ensayo debe ser >0 y <=100 mm.")
    surface = sample_patient(equipment.patient)
    base = calculate(events, equipment, stage, surface=surface).summary
    reference = np.asarray(base["psd_point_patient_mm"], dtype=float)
    original = np.asarray(equipment.patient.get("origin_on_table_mm", [0, 0, 0]), dtype=float)
    if original.shape != (3,) or not np.all(np.isfinite(original)):
        raise InputError("Posición del paciente sobre la mesa inválida.")
    scenarios = []
    for axis, label in enumerate(("x lateral", "y craneal", "z anterior")):
        for sign in (-1, 1):
            variation = np.zeros(3, dtype=float)
            variation[axis] = sign*offset_mm
            cfg = deepcopy(equipment.raw)
            cfg["patient"]["origin_on_table_mm"] = (original+variation).tolist()
            try:
                result = calculate(events, replace(equipment, raw=cfg), stage, surface=surface).summary
            except InputError as exc:
                if "Ningún rayo incide en la superficie" not in str(exc):
                    raise
                scenarios.append({"axis": label, "offset_mm": variation.tolist(),
                                  "status": "no_intersection", "psd_point_gy": None})
                continue
            scenarios.append({"axis": label, "offset_mm": variation.tolist(),
                              "status": "calculated", "psd_point_gy": result["psd_point_gy"],
                              "peak_position_patient_mm": result["psd_point_patient_mm"],
                              "peak_displacement_patient_mm": float(np.linalg.norm(
                                  np.asarray(result["psd_point_patient_mm"])-reference)),
                              "neighbourhood_1cm2_gy": result["psd_neighbourhood_1cm2"]["dose_gy"]
                              if result["psd_neighbourhood_1cm2"] else None})
    peaks = [case["psd_point_gy"] for case in scenarios if case["status"] == "calculated"]
    return {"stage": stage, "tested_offset_mm": offset_mm,
            "base_psd_point_gy": base["psd_point_gy"],
            "base_peak_position_patient_mm": reference.tolist(),
            "scenarios": scenarios,
            "psd_range_valid_scenarios_gy": [min(peaks), max(peaks)] if peaks else None,
            "relative_peak_change_range_valid_scenarios": [
                min(peaks)/base["psd_point_gy"]-1., max(peaks)/base["psd_point_gy"]-1.]
                if peaks else None,
            "interpretation": "Seis cambios hipotéticos de posición; el rango no es intervalo de confianza ni validación clínica.",
            "clinical_validation": False}
