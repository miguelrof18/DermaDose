"""Standalone, offline 3D map of the existing numerical dose grid.

The clipping controls only change the displayed triangles; they never change
the computed dose, the RDSR event data, or the mesh resolution.
"""
from __future__ import annotations

from importlib.resources import files
import json
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.io import to_html

from .engine import Calculation


def _sparse(values: np.ndarray) -> dict:
    rounded = np.round(np.asarray(values, dtype=float), 11)
    index = np.flatnonzero(rounded)
    return {"index": index.tolist(), "value": rounded[index].tolist()}


def _viewer_factor_audit(events: list[dict]) -> list[dict]:
    """Numeric node values per sampled beam, as read by the viewer panel."""
    audits = []
    for event in events:
        for audit in event.get("surface_factor_audit_per_sample") or []:
            if not audit:
                continue
            nodes = []
            for node in audit.get("source_nodes", []):
                try:
                    nodes.append({"k_bs": float(node["k_bs"]), "k_f": float(node["k_f"])})
                except (KeyError, TypeError, ValueError):
                    continue
            audits.append({"status": audit.get("status"), "source_nodes": nodes})
    return audits


def save_interactive(calculation: Calculation, path: Path) -> None:
    surface = calculation.surface
    is_mesh = surface.mesh is not None
    coords = np.asarray(surface.mesh.vertices if is_mesh else surface.points_mm, dtype=float)
    faces = np.asarray(surface.mesh.faces, dtype=int) if is_mesh else None
    total = np.asarray(calculation.total_gy, dtype=float).ravel()
    maximum = float(np.max(total))
    # Event maps are processed one at a time: peak and sparse encoding only,
    # so no events-by-surface matrix is held in memory for the viewer.
    sparse_modes, peaks = [_sparse(total)], [surface.points_mm[int(np.argmax(total))].tolist()]
    finite = bool(np.all(np.isfinite(total)))
    for event_map in calculation.event_maps_gy:
        event_map = np.asarray(event_map, dtype=float).ravel()
        finite &= bool(np.all(np.isfinite(event_map)))
        sparse_modes.append(_sparse(event_map))
        peaks.append(surface.points_mm[int(np.argmax(event_map))].tolist())
    if maximum <= 0 or not np.all(np.isfinite(coords)) or not finite:
        raise ValueError("Visor 3D: coordenadas o dosis no finitas.")
    # Dose-zero remains gray. Low positive dose starts in green, then moves
    # through yellow/orange/red to make small-but-real dose easier to see.
    palette_stops = [
        (0.00, "#2ca25f"),
        (0.33, "#ffe34d"),
        (0.67, "#fb8c00"),
        (1.00, "#d7191c"),
    ]

    def hex_to_rgb(color: str) -> np.ndarray:
        return np.asarray([int(color[i:i + 2], 16) for i in (1, 3, 5)], dtype=float)

    def rgb_to_hex(rgb: np.ndarray) -> str:
        values = np.clip(np.rint(rgb), 0, 255).astype(int)
        return "#{:02x}{:02x}{:02x}".format(*values)

    palette = []
    for i in range(256):
        fraction = i / 255
        for (x0, c0), (x1, c1) in zip(palette_stops, palette_stops[1:]):
            if fraction <= x1:
                local = 0 if x1 == x0 else (fraction - x0) / (x1 - x0)
                color = hex_to_rgb(c0) * (1 - local) + hex_to_rgb(c1) * local
                palette.append(rgb_to_hex(color))
                break
    no_dose = "#505d70"

    def colors(dose: np.ndarray) -> list[str]:
        levels = np.clip((dose/maximum*255).astype(int), 0, 255)
        return [palette[int(level)] if value > 0 else no_dose
                for level, value in zip(levels, dose)]

    first_color = colors(total)
    if is_mesh:
        fig = go.Figure(go.Mesh3d(
            x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
            i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
            facecolor=first_color, flatshading=False, hoverinfo="skip",
            lighting=dict(ambient=1.0, diffuse=0.0, specular=0.0, roughness=1.0),
            lightposition=dict(x=600, y=800, z=-1200), name="piel"))
    else:
        fig = go.Figure(go.Scatter3d(
            x=coords[:, 0], y=coords[:, 1], z=coords[:, 2], mode="markers",
            marker=dict(size=2.7, color=first_color, opacity=.85),
            hoverinfo="skip", name="puntos de piel"))
    fig.add_trace(go.Scatter3d(x=[peaks[0][0]], y=[peaks[0][1]], z=[peaks[0][2]],
                               mode="markers", marker=dict(size=7, color="#f9fafb",
                                                            line=dict(color="#151b29", width=2)),
                               hovertemplate="Máximo puntual: %{text}<extra></extra>",
                               text=[f"{maximum*1000:.3f} mGy"], name="máximo"))
    lower, upper = np.min(coords, axis=0), np.max(coords, axis=0)
    # Plane traces are initially hidden; JavaScript positions them at the sliders.
    for plane_color in ("#6bb8ff", "#f4ae6f", "#78e2b5"):
        fig.add_trace(go.Surface(x=[[0, 0], [0, 0]], y=[[0, 0], [0, 0]],
                                 z=[[0, 0], [0, 0]], showscale=False, opacity=.21,
                                 colorscale=[[0, plane_color], [1, plane_color]],
                                 hoverinfo="skip", visible=False, name="plano de corte"))
    # Stable indices used by the offline viewer: 5 rays, 6 source/isocentre,
    # 7 aperture, 8 inspectable centroids, 9 selection, 10 normal, 11 hotspots,
    # 12 schematic C-arm trajectory. Empty until JavaScript selects an event.
    for name, color, mode in (("rayos del haz", "#6ee7e0", "lines"),
                              ("fuente e isocentro", "#6ee7e0", "markers+text"),
                              ("campo de referencia", "#6ee7e0", "lines"),
                              ("inspección", "#ecf4ff", "markers"),
                              ("punto elegido", "#52d7ff", "markers"),
                              ("normal local", "#52d7ff", "lines"),
                              ("regiones calientes", "#ffcc70", "markers"),
                              ("arco esquemático", "#70849c", "lines")):
        fig.add_trace(go.Scatter3d(x=[], y=[], z=[], mode=mode, name=name,
                                  marker=dict(size=4, color=color), line=dict(width=3, color=color),
                                  textposition="top center", hoverinfo="skip", visible=False))
    ct_lps = calculation.summary.get("patient_coordinate_system") == "CT_LPS_mm"
    anchor = calculation.summary.get("xcat_anchor_surface_mm")
    if anchor is not None:
        fig.add_trace(go.Scatter3d(x=[anchor[0]], y=[anchor[1]], z=[anchor[2]], mode="markers",
                                   marker=dict(size=6, color="#00d5ff", symbol="diamond"),
                                   hovertemplate="Anclaje posterior de hombros XCAT<extra></extra>",
                                   name="Anclaje de hombros"))
        # The registered surface keeps LPS axes. The provisional surface uses
        # the engine axes, where the head points towards positive Y.
        longitudinal_axis = 2 if ct_lps else 1
        feet_label = np.asarray(anchor, dtype=float).copy()
        head_label = feet_label.copy()
        feet_label[longitudinal_axis] = lower[longitudinal_axis]
        head_label[longitudinal_axis] = upper[longitudinal_axis]
        labels = np.stack([feet_label, head_label])
        fig.add_trace(go.Scatter3d(
            x=labels[:, 0], y=labels[:, 1], z=labels[:, 2],
            mode="text", text=["Pies", "Cabeza · cabecera"],
            textfont=dict(color="#00d5ff", size=14), hoverinfo="skip",
            name="Orientación del XCAT",
        ))
    camera = (dict(eye=dict(x=.35, y=2.25, z=.20), up=dict(x=0, y=0, z=1)) if ct_lps else
              dict(eye=dict(x=.35, y=.20, z=-2.25), up=dict(x=0, y=1, z=0)))
    fig.update_layout(
        paper_bgcolor="#111722", font=dict(color="#e7eaf2", family="Arial, sans-serif"),
        margin=dict(l=0, r=0, t=0, b=0), showlegend=False,
        scene=dict(xaxis=dict(title="x lateral [mm]", backgroundcolor="#182131", color="#c6d0df"),
                   yaxis=dict(title="y posterior TC [mm]" if ct_lps else "y craneal [mm]", backgroundcolor="#182131", color="#c6d0df"),
                   zaxis=dict(title="z craneal TC [mm]" if ct_lps else "z anterior [mm]", backgroundcolor="#182131", color="#c6d0df"),
                   aspectmode="data", dragmode="orbit",
                   camera=camera),
        uirevision="dose-map")
    figure = to_html(fig, full_html=False, include_plotlyjs=True, div_id="dose-map",
                     config=dict(displaylogo=False, responsive=True,
                                 modeBarButtonsToRemove=["lasso2d", "select2d"]))
    payload = {
        "kind": surface.kind,
        "ct_lps": ct_lps,
        "vertices": np.round(coords, 3).tolist(),
        "faces": faces.tolist() if is_mesh else None,
        # Sparse per-mode encoding (non-zero faces only); the viewer expands it
        # to dense arrays on load. Values and rounding are unchanged.
        "sparse_doses": sparse_modes,
        "modes": ["Acumulada"] + [f"Evento {i + 1} ({r['plane']})"
                                  for i, r in enumerate(calculation.summary["events"])],
        "peaks": peaks, "peak_total_gy": maximum,
        "min_mm": lower.tolist(), "max_mm": upper.tolist(),
        "palette": palette, "no_dose_color": no_dose,
        "quantity": calculation.summary["quantity"],
        "centers_mm": np.round(surface.points_mm, 4).tolist(),
        "normals": np.round(surface.normals, 8).tolist(),
        "area_mm2": np.round(surface.area_mm2, 8).tolist(),
        "geometry": [report.get("geometry_samples", []) for report in calculation.summary["events"]],
        "metrics": calculation.summary.get("surface_metrics"),
        "factor_metadata": calculation.summary.get("factor_metadata", {}),
        "surface_factor_dataset": calculation.summary.get("surface_factor_dataset"),
        "surface_factor_audit": _viewer_factor_audit(calculation.summary.get("events", [])),
        "synthetic_reference": {"k_bs": 1.25, "k_f": 1.05, "identifier": "DEMO-BSF/DEMO-FSKIN"},
    }
    template = files("skindose_verifier").joinpath("data/visor_3d_template.html").read_text(encoding="utf-8")
    # JSON contains only finite numerical arrays and fixed text labels; escaping
    # angle brackets also avoids terminating the application/json script element.
    encoded = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).replace("<", "\\u003c")
    path.write_text(template.replace("__PLOT_HTML__", figure).replace("__VIEWER_DATA__", encoded),
                    encoding="utf-8")
