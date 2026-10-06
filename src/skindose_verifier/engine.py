"""Cálculo por evento y suma sobre coordenadas de paciente invariables."""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Sequence
import math
import time
from collections import OrderedDict

import numpy as np

from .config import Equipment
from .domain import InputError, IrradiationEvent
from .geometry import event_slices, project, beam_view, Projection
from .irp import verify_irp
from .patient import SkinSurface, sample_patient
from .metrics import neighbourhood_psd, dose_area_summary, surface_metrics
from .mesh_processing import mesh_quality


@dataclass
class Calculation:
    surface: SkinSurface
    total_gy: np.ndarray
    event_maps_gy: list[np.ndarray]
    summary: dict


def _factor_metadata(equipment: Equipment) -> dict:
    def describe(factor) -> dict:
        data = {
            "identifier": factor.identifier,
            "provenance": factor.provenance,
            "relative_standard_uncertainty": factor.relative_standard_uncertainty,
            "representation": "constant" if factor.value is not None else "tabulated_interpolated",
        }
        if factor.value is not None:
            data["value"] = factor.value
        if factor.axes is not None:
            data["axes"] = {name: {"minimum": float(min(values)), "maximum": float(max(values)),
                                   "points": [float(v) for v in values]}
                            for name, values in factor.axes.items()}
        return data
    metadata = {name: describe(factor) for name, factor in equipment.factors.items()}
    by_mode = {}
    for mode, factors in equipment.factors_by_mode.items():
        by_mode[mode] = {name: describe(factor) for name, factor in factors.items()}
    return {"default": metadata, "by_mode": by_mode}


def normalize(ev: IrradiationEvent, equipment: Equipment) -> IrradiationEvent:
    """Only explicitly versioned device geometry may fill missing distances."""
    g = equipment.geometry_for(ev.plane)
    updates = {}
    q = ev.quality.copy()
    for field, cfg in (("sid_mm", "source_to_detector_mm"), ("sod_mm", "source_to_isocenter_mm"),
                       ("srp_mm", "source_to_irp_mm")):
        if getattr(ev, field) is None:
            updates[field] = float(g[cfg])
            q[field] = "defaulted_from_equipment_configuration"
    system = g["table_coordinates"]
    if ev.table_system != system:
        keys = ("table_x_iso", "table_y_iso", "table_z_iso") if system == "isocenter" else (
            "table_lat", "table_long", "table_height")
        raw = ev.original
        if all(k in raw for k in keys):
            updates["table_mm"] = tuple(raw[k][0]["value"] for k in keys)
            updates["table_system"] = system
            q["table_mm"] = f"selected_{system}_from_rdsr"
    updates["quality"] = q
    return replace(ev, **updates)


def _factor(factors: dict, name: str, context: dict) -> np.ndarray:
    return factors[name].evaluate(**context)


def _dose(ev: IrradiationEvent, equipment: Equipment, surface: SkinSurface,
          stage: str, geo: Projection | None = None) -> tuple[np.ndarray, dict]:
    geo = project(ev, equipment, surface) if geo is None else geo
    count = len(surface.points_mm)
    dose = np.zeros(count)
    if not np.any(geo.visible_mask):
        return dose, {"irradiated_cells": 0, "table_cells": 0, "maximum_gy": 0., "maximum_location_mm": None}
    mask = geo.visible_mask
    distance = geo.distance_mm[mask]
    components: dict[str, np.ndarray] = {}
    kerma = float(ev.ka_r_gy) * (geo.beam.irp_distance_mm/distance)**2
    components["inverse_square"] = (geo.beam.irp_distance_mm/distance)**2
    correction = np.ones(len(distance))
    if stage != "free_air":
        factors = equipment.factors_for(ev)
        context = {"kvp": ev.kvp, "filtration_mm_cu": ev.filtration_mm_cu,
                   "field_area_cm2": geo.field_area_cm2[mask],
                   "table_angle_deg": geo.table_angle_deg[mask],
                   "skin_angle_deg": geo.skin_angle_deg[mask],
                   "field_fraction_x": np.clip(geo.field_fraction_xy[mask, 0], -1., 1.),
                   "field_fraction_y": np.clip(geo.field_fraction_xy[mask, 1], -1., 1.),
                   "patient_thickness_mm": equipment.patient.get("thickness_mm",
                       2*float(equipment.patient.get("semi_axis_z_mm", 0.)))}
        cf = _factor(factors, "cf", context)
        correction *= cf
        components["cf"] = cf
        # Optional measured profile, relative to central-axis kerma; never a
        # guessed heel correction. See TG357 and Gañán Mora et al. (2024).
        if "beam_profile" in factors:
            profile = _factor(factors, "beam_profile", context)
            correction *= profile
            components["beam_profile"] = profile
        table = geo.table_intersection_mask[mask]
        if np.any(table):
            table_context = {key: value[table] if isinstance(value, np.ndarray) else value
                             for key, value in context.items()}
            if "table_total" in factors:
                transmission = _factor(factors, "table_total", table_context)
            else:
                transmission = _factor(factors, "table_normal", table_context)*_factor(
                    factors, "table_angular", table_context)
            correction[table] *= transmission
            components["table_transmission"] = transmission
        if stage == "skin_dose":
            bsf = _factor(factors, "bsf", context)
            fskin = _factor(factors, "f_skin", context)
            correction *= bsf*fskin
            components["bsf"] = bsf
            components["f_skin"] = fskin
    dose[mask] = kerma * correction
    idx = int(np.argmax(dose))
    return dose, {
        "irradiated_cells": int(np.sum(mask)),
        "irradiated_area_cm2": float(np.sum(surface.area_mm2[mask])/100.),
        "table_cells": int(np.sum(geo.table_intersection_mask)),
        "maximum_gy": float(dose[idx]),
        "maximum_location_mm": surface.points_mm[idx].tolist(),
        "component_ranges": {name: [float(np.min(values)), float(np.max(values))]
                             for name, values in components.items()},
    }


def uncertainty(equipment: Equipment, stage: str, events: Sequence[IrradiationEvent]) -> dict:
    names = [] if stage == "free_air" else ["cf"] + sorted({
        name for ev in events for name in equipment.factors_for(ev)
        if name in ("table_total", "table_normal", "table_angular", "beam_profile")})
    if stage == "skin_dose":
        names += ["bsf", "f_skin"]
    known: dict[str, float] = {}
    missing = ["posición real del paciente", "geometría real del campo", "incertidumbre de discretización"]
    for name in names:
        uncertainties = [equipment.factors_for(ev)[name].relative_standard_uncertainty
                         for ev in events if name in equipment.factors_for(ev)]
        if any(u is None for u in uncertainties):
            missing.append(f"incertidumbre de {name}")
        else:
            known[name] = max(uncertainties)
    u_kar = equipment.raw.get("rdsr_relative_standard_uncertainty")
    if u_kar is None:
        missing.append("incertidumbre de Ka,r informado")
    else:
        known["ka_r"] = float(u_kar)
    return {"partial_relative_standard_uncertainty": math.sqrt(sum(x*x for x in known.values())),
            "quantified_components": known, "unquantified_components": missing,
            "complete_budget": False}


def calculate(events: Sequence[IrradiationEvent], equipment: Equipment, stage: str = "free_air", *,
              surface: SkinSurface | None = None, detailed: bool = True,
              cache_geometry: bool = True) -> Calculation:
    """Evaluate a fixed surface. Reuse `surface` only for unchanged local anatomy.

    Bounded, call-local cache cannot leak geometry between configurations.
    `detailed=False` omits costly local-area/cluster reports for parameter draws.
    """
    started = time.perf_counter()
    if not events:
        raise InputError("No hay eventos de irradiación.")
    surface = sample_patient(equipment.patient) if surface is None else surface
    surface_seconds = time.perf_counter()-started
    cache: OrderedDict[tuple, Projection] = OrderedDict()
    cache_hits = 0
    total_rays = total_pairs = 0
    maps = []
    reports = []
    for raw in events:
        event_started = time.perf_counter()
        ev = normalize(raw, equipment)
        equipment.validate_event(ev, stage)
        sliced = event_slices(ev, equipment)
        event_map = np.zeros(len(surface.points_mm))
        subreports = []
        event_rays = event_pairs = 0
        for piece in sliced:
            key = (piece.plane, piece.primary_deg, piece.secondary_deg, tuple(piece.table_mm),
                   piece.sid_mm, piece.sod_mm, piece.srp_mm, piece.field_width_mm, piece.field_height_mm)
            geo = cache.get(key) if cache_geometry else None
            if geo is None:
                geo = project(piece, equipment, surface)
                event_rays += geo.statistics["rays_traced"]
                event_pairs += geo.statistics["candidate_triangle_pairs"]
                if cache_geometry:
                    cache[key] = geo
                    if len(cache) > 2:
                        cache.popitem(last=False)
            else:
                cache_hits += 1
                cache.move_to_end(key)
            piece_map, sub = _dose(piece, equipment, surface, stage, geo)
            event_map += piece_map
            subreports.append(sub)
        idx = int(np.argmax(event_map))
        illuminated = bool(np.any(event_map))
        total_rays += event_rays
        total_pairs += event_pairs
        reports.append({"event_uid": ev.uid, "plane": ev.plane, "ka_r_gy": ev.ka_r_gy,
                        "kap_gy_m2": ev.kap_gy_m2, "status": ev.quality,
                        "sample_count": len(sliced),
                        "geometry_samples": [beam_view(piece, equipment) for piece in sliced],
                        "performance": {"seconds": time.perf_counter()-event_started,
                                        "rays_traced": event_rays, "candidate_triangle_pairs": event_pairs},
                        "uniform_sweep_assumed": len(sliced)>1,
                        "maximum_gy": float(event_map[idx]),
                        "maximum_patient_point_mm": surface.points_mm[idx].tolist() if illuminated else None,
                        "irradiated_cells": int(np.count_nonzero(event_map)),
                        "irradiated_area_cm2_approx": float(np.sum(surface.area_mm2[event_map > 0])/100.),
                        "no_skin_intersection_warning": not illuminated,
                        "table_intersection_cells_per_sample": [x["table_cells"] for x in subreports],
                        "irp_verification": verify_irp(ev, equipment),
                        "factor_versions": {name: f.identifier for name, f in equipment.factors_for(ev).items()}
                        if stage != "free_air" else {},
                        "factor_ranges_per_sample": [x.get("component_ranges", {}) for x in subreports]})
        maps.append(event_map)
    total = np.sum(np.stack(maps), axis=0)
    if not np.any(total > 0):
        raise InputError("Ningún rayo incide en la superficie de piel modelada; revisar campo, geometría y paciente.")
    peak_index = int(np.argmax(total))
    peak_event_contributions = []
    for report, event_map in zip(reports, maps):
        value = float(event_map[peak_index])
        peak_event_contributions.append({
            "event_uid": report["event_uid"],
            "dose_at_global_psd_point_gy": value,
            "fraction_of_global_psd": value/float(total[peak_index]) if total[peak_index] > 0 else 0.,
            "event_own_maximum_gy": report["maximum_gy"],
        })
    peak_event_contributions.sort(key=lambda item: item["dose_at_global_psd_point_gy"], reverse=True)
    for rank, item in enumerate(peak_event_contributions, start=1):
        item["rank"] = rank
    summary = {
        "status": "research_unvalidated_no_clinical_use",
        "quantity": {"free_air": "free_in_air_kerma_at_skin_Gy",
                     "calibrated_air": "calibrated_incident_air_kerma_at_skin_Gy",
                     "skin_dose": f"entrance_absorbed_dose_in_{equipment.raw.get('dose_medium')}_Gy"}[stage],
        "stage": stage, "events": reports, "event_count": len(reports),
        "total_rdsr_ka_r_gy": float(sum(ev.ka_r_gy for ev in events)),
        "psd_point_gy": float(total[peak_index]),
        "psd_point_patient_mm": surface.points_mm[peak_index].tolist(),
        "global_psd_event_contributions": peak_event_contributions,
        "psd_neighbourhood_1cm2": neighbourhood_psd(surface, total) if detailed else None,
        "surface_dose_area_approx": dose_area_summary(surface, total),
        "spatial_resolution_target_mm": equipment.patient["resolution_mm"] if surface.kind != "mesh" else None,
        "surface_metrics": surface_metrics(surface, total, equipment.raw.get("dose_area_thresholds_gy")) if detailed else None,
        "mesh_geometry": mesh_quality(surface.mesh) if surface.mesh is not None and detailed else None,
        "uncertainty": uncertainty(equipment, stage, [normalize(ev, equipment) for ev in events]),
        "events_outside_model": [r["event_uid"] for r in reports if r["no_skin_intersection_warning"]],
        "calibration_versions": {k: f.identifier for k, f in equipment.factors.items()},
        "factor_metadata": _factor_metadata(equipment),
        "configuration_version": equipment.version,
        "assumptions": ["Colimación rectangular; perfil uniforme salvo beam_profile calibrado", "Paciente rígido y centrado según configuración",
                        "Mapa calculado sobre puntos de superficie; sin transporte Monte Carlo"] + (
                            ["Malla de referencia, no anatomía individual; resolución desigual de caras"]
                            if surface.kind == "mesh" else []),
    }
    summary["performance"] = {"total_seconds": time.perf_counter()-started, "surface_preparation_seconds": surface_seconds,
                              "projection_cache_hits": cache_hits, "rays_traced": total_rays,
                              "candidate_triangle_pairs": total_pairs, "surface_samples": len(total),
                              "timing_excludes_export": True}
    return Calculation(surface, total, maps, summary)
