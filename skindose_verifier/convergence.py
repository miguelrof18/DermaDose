"""Compare spatial resolutions without treating agreement as clinical validation."""
from copy import deepcopy
from dataclasses import replace
import math
import time
import numpy as np

from .config import Equipment
from .domain import InputError, IrradiationEvent
from .engine import calculate
from .patient import sample_patient


def compare_resolutions(events: list[IrradiationEvent] | tuple[IrradiationEvent, ...],
                        equipment: Equipment, stage: str, resolutions_mm: list[float]) -> dict:
    if len(resolutions_mm) < 2 or any(not math.isfinite(r) or r <= 0 for r in resolutions_mm):
        raise InputError("Especificar al menos dos resoluciones positivas en mm.")
    if len(set(resolutions_mm)) != len(resolutions_mm):
        raise InputError("Resoluciones repetidas.")
    values = []
    is_mesh = equipment.patient["type"] == "mesh"
    for step in sorted(resolutions_mm, reverse=True):
        cfg = deepcopy(equipment.raw)
        cfg["patient"]["mesh_max_edge_mm" if is_mesh else "resolution_mm"] = step
        started = time.perf_counter()
        calc = calculate(events, replace(equipment, raw=cfg), stage, detailed=not is_mesh)
        result = calc.summary
        values.append({"resolution_mm": step, "psd_point_gy": result["psd_point_gy"],
                       "point_mm": result["psd_point_patient_mm"],
                       "psd_neighbourhood_1cm2": result["psd_neighbourhood_1cm2"],
                       "surface_samples": len(calc.total_gy),
                       "irradiated_area_cm2": result["surface_dose_area_approx"]["irradiated_area_cm2"],
                       "dose_area_integral_gy_cm2": float(np.dot(calc.total_gy, calc.surface.area_mm2)/100),
                       "seconds": time.perf_counter()-started,
                       "maximum_triangle_area_mm2": float(calc.surface.area_mm2.max()) if is_mesh else None})
    coarse, fine = values[-2:]
    displacement = math.dist(coarse["point_mm"], fine["point_mm"])
    return {"stage": stage, "values": values,
            "method": "shared_edge_refinement_fixed_anatomy" if is_mesh else "resampled_analytic_surface",
            "mesh_refinement_roi_mm": equipment.patient.get("mesh_refinement_roi_mm") if is_mesh else None,
            "relative_integral_change_last_refinement": fine["dose_area_integral_gy_cm2"]/coarse["dose_area_integral_gy_cm2"]-1,
            "relative_irradiated_area_change_last_refinement": fine["irradiated_area_cm2"]/coarse["irradiated_area_cm2"]-1,
            "relative_point_psd_change_last_refinement":
                (fine["psd_point_gy"]-coarse["psd_point_gy"])/coarse["psd_point_gy"]
                if coarse["psd_point_gy"] else None,
            "point_displacement_mm_last_refinement": displacement,
            "clinical_validation": False,
            "interpretation": "Sensibilidad de discretización; no demuestra exactitud anatómica ni convergencia monótona en bordes de campo."}
