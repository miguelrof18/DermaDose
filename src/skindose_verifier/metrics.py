"""Surface-dose metrics: area weighting and adjacency, no volumetric DVH.

The neighbourhood is a discrete approximation, not an exact 1 cm² dosimeter.
"""
from __future__ import annotations
import heapq
import math
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from .domain import InputError
from .patient import SkinSurface


def neighbourhood_psd(surface: SkinSurface, dose: np.ndarray, area_mm2: float = 100.) -> dict | None:
    """Area-weighted local disc in surface coordinates; repeated theta for ellipse."""
    if area_mm2 <= 0:
        raise InputError("Área de promediado no positiva.")
    if surface.kind == "mesh":
        return _mesh_neighbourhood_psd(surface, dose, area_mm2)
    if np.max(surface.area_mm2) > area_mm2:
        return None  # spatial resolution is too coarse to resolve requested averaging area
    y = surface.axes[0]
    if surface.kind == "plane":
        x = surface.axes[1]
        xx, yy = np.meshgrid(x, y)
        uv = np.column_stack((xx.ravel(), yy.ravel()))
        periodic = False
    else:
        theta = surface.axes[1]
        n = len(theta)
        # Arc-length coordinate from the discrete representative areas.
        ds = surface.area_mm2[:n]/(y[1]-y[0] if len(y)>1 else surface.area_mm2[0]/1.)
        arc = np.cumsum(ds) - ds/2
        circumference = float(np.sum(ds))
        ss, yy = np.meshgrid(arc, y)
        uv = np.column_stack((ss.ravel(), yy.ravel()))
        periodic = True
    if periodic:
        tiled = np.concatenate((uv, uv+[-circumference, 0.], uv+[circumference, 0.]))
        tree = cKDTree(tiled)
    else:
        tree = cKDTree(uv)
    radius = math.sqrt(area_mm2/math.pi)
    result = np.full(len(dose), -np.inf)
    for i, ids in enumerate(tree.query_ball_point(uv, radius)):
        ids = np.unique(np.asarray(ids, dtype=int) % len(dose))
        areas = surface.area_mm2[ids]
        # Boundary discs missing >15% of intended surface do not define a PSD.
        if .85*area_mm2 <= areas.sum() <= 1.25*area_mm2:
            result[i] = float(np.dot(dose[ids], areas)/areas.sum())
    if not np.any(np.isfinite(result)):
        return None
    idx = int(np.argmax(result))
    selected = np.unique(np.asarray(tree.query_ball_point(uv[idx], radius), dtype=int) % len(dose))
    return {"dose_gy": float(result[idx]), "patient_point_mm": surface.points_mm[idx].tolist(),
            "sampled_area_mm2": float(surface.area_mm2[selected].sum()),
            "target_area_mm2": area_mm2,
            "definition": "media ponderada en vecindad circular discreta; no equivale necesariamente a PSD clínica"}


def _mesh_neighbourhood_psd(surface: SkinSurface, dose: np.ndarray, area_mm2: float) -> dict | None:
    """Local neighbourhood along adjacent triangles; never jump between limbs."""
    assert surface.mesh is not None
    radius = math.sqrt(area_mm2/math.pi)
    graph = surface.face_graph
    best_value, best_idx, best_area = -1., -1, 0.
    for origin in np.flatnonzero(dose > 0):
        origin = int(origin)
        distances = {origin: 0.}
        queue = [(0., origin)]
        indices = []
        while queue:
            length, node = heapq.heappop(queue)
            if length > radius:
                break
            if length > distances[node]:
                continue
            indices.append(node)
            start, stop = graph.indptr[node:node+2]
            for neighbour, segment in zip(graph.indices[start:stop], graph.data[start:stop]):
                candidate = length + segment
                if candidate <= radius and candidate < distances.get(neighbour, float("inf")):
                    distances[neighbour] = candidate
                    heapq.heappush(queue, (candidate, neighbour))
        weights = surface.area_mm2[indices]
        area = float(weights.sum())
        if not (.85*area_mm2 <= area <= 1.25*area_mm2):
            continue
        mean = float(np.dot(dose[indices], weights)/area)
        if mean > best_value:
            best_value, best_idx, best_area = mean, origin, area
    if best_idx < 0:
        return None
    return {"dose_gy": best_value, "patient_point_mm": surface.points_mm[best_idx].tolist(),
            "sampled_area_mm2": best_area, "target_area_mm2": area_mm2,
            "definition": "media ponderada por cara en disco geodésico aproximado sobre caras adyacentes"}


def dose_area_summary(surface: SkinSurface, dose: np.ndarray) -> dict:
    """Descriptive surface areas from sampled cells, never an organ dose or DVH."""
    peak = float(np.max(dose))
    areas = surface.area_mm2
    return {
        "unit": "cm2", "sampling": "suma de áreas de celdas según dosis en su punto representativo",
        "irradiated_area_cm2": float(np.sum(areas[dose > 0])/100.),
        "area_at_least_10_percent_peak_cm2": float(np.sum(areas[(dose > 0) & (dose >= peak*.10)])/100.),
        "area_at_least_50_percent_peak_cm2": float(np.sum(areas[(dose > 0) & (dose >= peak*.50)])/100.),
        "area_at_least_90_percent_peak_cm2": float(np.sum(areas[(dose > 0) & (dose >= peak*.90)])/100.),
    }


def weighted_dose_metrics(dose: np.ndarray, area_mm2: np.ndarray) -> dict:
    """D_p is the dose exceeded by at least p% of the declared area.

    Uses a step-function inverse of the area CDF; does not interpolate across
    zero-dose skin or give every triangle the same statistical weight.
    """
    if len(dose) == 0:
        return {"area_cm2": 0., "dmean_gy": None, "dmax_gy": None,
                "d2_percent_gy": None, "d5_percent_gy": None, "d10_percent_gy": None}
    order = np.argsort(-dose, kind="stable")
    cum = np.cumsum(area_mm2[order])/area_mm2.sum()
    out = {"area_cm2": float(area_mm2.sum()/100), "dmean_gy": float(np.average(dose, weights=area_mm2)),
           "dmax_gy": float(dose.max())}
    for p in (2, 5, 10):
        index = min(len(dose)-1, int(np.searchsorted(cum, p/100, side="left")))
        out[f"d{p}_percent_gy"] = float(dose[order[index]])
    return out


def hotspot_regions(surface: SkinSurface, dose: np.ndarray, fraction: float = .9) -> dict:
    if not 0 < fraction <= 1:
        raise InputError("Fracción de hotspot debe estar en (0,1].")
    threshold = float(dose.max()*fraction)
    ids = np.flatnonzero((dose > 0) & (dose >= threshold))
    if surface.mesh is None:
        return {"threshold_gy": threshold, "available": False,
                "reason": "Clustering por adyacencia implementado para mallas triangulares."}
    n, labels = connected_components(surface.face_graph[ids][:, ids], directed=False)
    regions = []
    for k in range(n):
        region = ids[labels == k]
        peak = int(region[np.argmax(dose[region])])
        areas = surface.area_mm2[region]
        regions.append({"face_count": len(region), "area_cm2": float(areas.sum()/100),
                        "dmean_gy": float(np.average(dose[region], weights=areas)),
                        "dmax_gy": float(dose[peak]), "peak_face_index": peak,
                        "peak_patient_mm": surface.points_mm[peak].tolist(),
                        "area_centroid_patient_mm": np.average(surface.points_mm[region], axis=0, weights=areas).tolist()})
    regions.sort(key=lambda item: (-item["dmax_gy"], -item["area_cm2"]))
    return {"available": True, "threshold_gy": threshold, "fraction_of_peak": fraction,
            "connectivity": "shared_triangle_edge", "region_count": n,
            "regions": regions[:100], "regions_truncated": n > 100}


def surface_metrics(surface: SkinSurface, dose: np.ndarray,
                    thresholds_gy: list[float] | None = None) -> dict:
    if dose.shape != surface.area_mm2.shape or not np.all(np.isfinite(dose)) or np.any(dose < 0):
        raise InputError("Mapa de dosis inválido para métricas.")
    thresholds = np.asarray(thresholds_gy if thresholds_gy is not None else [.001, .01, .02, .05, .1, .5, 1., 2., 3., 5.])
    if thresholds.ndim != 1 or not np.all(np.isfinite(thresholds)) or np.any(thresholds < 0):
        raise InputError("Umbrales de dosis no válidos.")
    positive = dose > 0
    maximum = float(dose.max())
    levels = np.linspace(0, maximum, 101) if maximum else np.array([0.])
    area = surface.area_mm2
    return {"definition": "DAH superficial; no DVH volumétrico ni predictor de toxicidad",
            "whole_model_surface": weighted_dose_metrics(dose, area),
            "irradiated_surface": weighted_dose_metrics(dose[positive], area[positive]),
            "area_above_thresholds": [{"threshold_gy": float(t), "area_cm2": float(area[dose >= t].sum()/100)}
                                      for t in sorted(set(thresholds))],
            "cumulative_dose_area_histogram": {"dose_gy": levels.tolist(),
                "area_at_least_cm2": [float(area[dose >= t].sum()/100) for t in levels]},
            "dose_area_integral_gy_cm2": float(np.dot(dose, area)/100),
            "integral_is_not_kap_or_energy": True,
            "hotspots": hotspot_regions(surface, dose)}

