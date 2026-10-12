"""Estimaciones independientes o controles de consistencia para Ka,r."""
from __future__ import annotations

from .config import Equipment, Factor
from .domain import IrradiationEvent, InputError


def verify_irp(ev: IrradiationEvent, equipment: Equipment) -> dict:
    """Tube yield is independent of reported Ka,r; KAP/area is only a crosscheck."""
    g = equipment.geometry_for(ev.plane)
    srp = float(ev.srp_mm)
    out: dict = {"reported_ka_r_gy": ev.ka_r_gy, "tube_output": None, "kap_area": None}
    model = equipment.raw.get("tube_output")
    if model and ev.exposure_mas is not None:
        dose_ref = float(model["reference_distance_mm"])
        if dose_ref <= 0:
            raise InputError("Distancia de referencia del output inválida.")
        yield_factor = Factor.from_dict("tube_yield_gy_per_mas", model["yield_gy_per_mas"])
        context = {"kvp": ev.kvp, "filtration_mm_cu": ev.filtration_mm_cu}
        missing = [name for name in (yield_factor.axes or {}) if context.get(name) is None]
        if missing:
            out["tube_output"] = {"classification": "unavailable_missing_yield_parameters",
                                  "missing": missing}
        else:
            value = float(yield_factor.evaluate(**context))
            independent = value * ev.exposure_mas * (dose_ref/srp)**2
            out["tube_output"] = {
                "classification": "independent_if_yield_measured_and_mas_credible",
                "ka_r_gy": independent, "yield_identifier": yield_factor.identifier,
                "relative_difference_to_reported": None if not ev.ka_r_gy else independent/ev.ka_r_gy-1.}
    if ev.kap_gy_m2 is not None and ev.field_width_mm is not None and ev.field_height_mm is not None:
        field_ref_distance = float(ev.sid_mm) if g["field_reference"] == "detector" else float(ev.sod_mm)
        area_m2 = float(ev.field_width_mm)*float(ev.field_height_mm)*(
            srp/field_ref_distance)**2/1e6
        if area_m2 <= 0:
            raise InputError("Área del campo en IRP no positiva.")
        estimate = ev.kap_gy_m2 / area_m2
        out["kap_area"] = {
            "classification": "consistency_only_assumes_uniform_beam_and_known_kap_plane",
            "ka_axis_gy": estimate, "area_at_irp_m2": area_m2,
            "relative_difference_to_reported": None if not ev.ka_r_gy else estimate/ev.ka_r_gy-1.}
    return out
