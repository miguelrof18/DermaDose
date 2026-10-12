"""Seeded parameter propagation; distribution assumptions are user inputs.

JCGM 101:2008 (doi:10.59161/JCGM101-2008) motivates propagation through
the complete model. Stratified draws use scipy.stats.qmc.LatinHypercube.
This is NOT photon transport Monte Carlo or a complete uncertainty budget.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, replace
import time
from typing import Sequence

import numpy as np
from scipy.stats import qmc, norm

from .config import Equipment
from .domain import InputError, IrradiationEvent
from .engine import calculate, normalize
from .patient import sample_patient


@dataclass
class PropagationResult:
    summary: dict
    points_mm: np.ndarray
    area_mm2: np.ndarray
    mean_gy: np.ndarray
    standard_deviation_gy: np.ndarray
    probability_exceeding_threshold: np.ndarray


PARAMETERS = {"patient_x_mm", "patient_y_mm", "patient_z_mm", "source_x_mm", "source_y_mm", "source_z_mm",
              "primary_deg", "secondary_deg", "sod_mm", "ka_scale", "cf_scale", "bsf_scale",
              "f_skin_scale", "table_scale", "profile_scale"}


def draw_parameters(specification: dict) -> tuple[list[str], np.ndarray]:
    """Independent dimensions, but each sampled bias is shared across events."""
    params = specification.get("parameters", {})
    n, seed = specification.get("samples", 64), specification.get("seed", 20260927)
    if not isinstance(n, int) or not 8 <= n <= 4096 or not isinstance(seed, int) or seed < 0:
        raise InputError("Incertidumbre: samples entero entre 8 y 4096; seed entero no negativo.")
    if not params or set(params)-PARAMETERS:
        raise InputError("Incertidumbre: parámetros ausentes o no admitidos.")
    if not specification.get("provenance"):
        raise InputError("Declarar provenance de las distribuciones, incluso si son hipotéticas.")
    method = specification.get("method", "lhs")
    if method not in ("lhs", "monte_carlo"):
        raise InputError("Método de propagación admitido: lhs o monte_carlo.")
    names = list(params)
    unit = (qmc.LatinHypercube(d=len(names), seed=seed).random(n) if method == "lhs"
            else np.random.default_rng(seed).random((n, len(names))))
    out = np.zeros_like(unit)
    for j, name in enumerate(names):
        p = params[name]
        mean = float(p.get("mean", 1. if name.endswith("_scale") else 0.))
        sd = float(p["standard_deviation"])
        if not np.isfinite(mean+sd) or sd < 0:
            raise InputError(f"Distribución inválida para {name}.")
        distribution = p.get("distribution")
        if distribution == "normal":
            out[:, j] = mean + sd*norm.ppf(np.clip(unit[:, j], 1e-12, 1-1e-12))
        elif distribution == "uniform":
            out[:, j] = mean + np.sqrt(3)*sd*(2*unit[:, j]-1)
        elif distribution == "lognormal" and name.endswith("_scale") and mean > 0:
            sigma = np.sqrt(np.log1p((sd/mean)**2))
            out[:, j] = np.exp(np.log(mean)-sigma**2/2 + sigma*norm.ppf(np.clip(unit[:, j], 1e-12, 1-1e-12)))
        else:
            raise InputError(f"Distribución no admitida para {name}.")
        if name.endswith("_scale") and np.any(out[:, j] <= 0):
            raise InputError(f"{name}: una realización no positiva; utilizar distribución físicamente válida.")
    return names, out


def _scaled_factor(data: dict, scale: float) -> None:
    if "value" in data:
        data["value"] *= scale
    else:
        data["grid"] = (np.asarray(data["grid"])*scale).tolist()


def propagate_parameters(events: Sequence[IrradiationEvent], equipment: Equipment,
                         stage: str, specification: dict) -> PropagationResult:
    started = time.perf_counter()
    names, draws = draw_parameters(specification)
    if stage == "free_air" and any(name in names for name in ("cf_scale", "bsf_scale", "f_skin_scale", "table_scale", "profile_scale")):
        raise InputError("Factores dosimétricos no se aplican en free_air.")
    if stage != "skin_dose" and any(name in names for name in ("bsf_scale", "f_skin_scale")):
        raise InputError("BSF/f_skin requieren etapa skin_dose.")
    surface = sample_patient(equipment.patient)
    evs = [normalize(e, equipment) for e in events]
    base = calculate(evs, equipment, stage, surface=surface, detailed=False)
    means, m2, exceed = (np.zeros(len(surface.points_mm)) for _ in range(3))
    threshold = float(specification.get("threshold_gy", .02))
    if not np.isfinite(threshold) or threshold <= 0:
        raise InputError("threshold_gy debe ser positivo y finito.")
    mapping = {"cf_scale": ("cf",), "bsf_scale": ("bsf",), "f_skin_scale": ("f_skin",),
               "table_scale": ("table_total", "table_normal"), "profile_scale": ("beam_profile",)}
    peaks, locations = [], []
    for i, values in enumerate(draws):
        params = dict(zip(names, values))
        cfg = deepcopy(equipment.raw)
        origin = np.asarray(cfg["patient"].get("origin_on_table_mm", [0, 0, 0]), dtype=float)
        cfg["patient"]["origin_on_table_mm"] = (origin + [params.get(f"patient_{axis}_mm", 0.) for axis in "xyz"]).tolist()
        # A source/isocenter frame translation: detector follows, preserving
        # collimation direction and source-reference distances. No independent
        # source motion relative to a fixed detector is implied.
        delta_source = np.array([params.get(f"source_{axis}_mm", 0.) for axis in "xyz"])
        if np.any(delta_source):
            g = cfg["geometry"]
            base_iso = np.asarray(g.get("isocenter_mm", [0, 0, 0]))
            for plane in g.get("planes", {}).values():
                plane["isocenter_mm"] = (np.asarray(plane.get("isocenter_mm", base_iso))+delta_source).tolist()
            g["isocenter_mm"] = (base_iso+delta_source).tolist()
        for name, factor_names in mapping.items():
            if name not in params:
                continue
            matched = False
            for block in [cfg["factors"], *cfg.get("factors_by_mode", {}).values()]:
                for factor_name in factor_names:
                    if factor_name in block:
                        _scaled_factor(block[factor_name], params[name]); matched = True
            if not matched:
                raise InputError(f"{name}: ningún factor correspondiente configurado.")
        changed = [replace(e.scaled(params.get("ka_scale", 1.)),
                           primary_deg=e.primary_deg+params.get("primary_deg", 0.),
                           secondary_deg=e.secondary_deg+params.get("secondary_deg", 0.),
                           primary_end_deg=None if e.primary_end_deg is None else e.primary_end_deg+params.get("primary_deg", 0.),
                           secondary_end_deg=None if e.secondary_end_deg is None else e.secondary_end_deg+params.get("secondary_deg", 0.),
                           sod_mm=e.sod_mm+params.get("sod_mm", 0.)) for e in evs]
        try:
            result = calculate(changed, Equipment.from_dict(cfg), stage, surface=surface, detailed=False)
        except InputError as exc:
            raise InputError(f"Propagación interrumpida en realización {i+1}; no se descartan muestras inválidas: {exc}") from exc
        dose = result.total_gy
        difference = dose-means
        means += difference/(i+1)
        m2 += difference*(dose-means)
        exceed += dose >= threshold
        peaks.append(result.summary["psd_point_gy"])
        locations.append(result.summary["psd_point_patient_mm"])
    peaks = np.asarray(peaks)
    summary = {"specification": specification, "parameter_order": names,
               "parameter_samples": draws.tolist(), "base_psd_gy": base.summary["psd_point_gy"],
               "psd_samples_gy": peaks.tolist(), "peak_positions_patient_mm": locations,
               "psd_mean_gy": float(peaks.mean()), "psd_sample_standard_deviation_gy": float(peaks.std(ddof=1)),
               "psd_quantiles_025_50_975_gy": np.quantile(peaks, [.025, .5, .975]).tolist(),
               "threshold_gy": threshold, "samples": len(draws), "seconds": time.perf_counter()-started,
               "dependence": "independent parameters; shared systematic bias across all events in each draw",
               "conditional_model_distribution": True, "complete_uncertainty_budget": False,
               "photon_transport_monte_carlo": False, "clinical_validation": False,
               "excluded": ["mesh discretization and anatomy/segmentation error", "unmodelled scatter/transport",
                            "calibration interpolation error", "temporal random error within each event",
                            "input parameter correlations"],
               "interpretation": "Cuantiles condicionados a las distribuciones introducidas; no son un intervalo clínico validado. Comprobar estabilidad con más muestras y semillas."}
    return PropagationResult(summary, surface.points_mm, surface.area_mm2, means,
                             np.sqrt(np.maximum(m2/(len(draws)-1), 0)), exceed/len(draws))
