"""Standalone, offline 3D map of the existing numerical dose grid.

The clipping controls only change the displayed triangles; they never change
the computed dose, the RDSR event data, or the mesh resolution.
"""
from __future__ import annotations

from importlib.resources import files
import json
from pathlib import Path

import matplotlib
import numpy as np
import plotly.graph_objects as go
from plotly.io import to_html

from .engine import Calculation


def save_interactive(calculation: Calculation, path: Path) -> None:
    surface = calculation.surface
    is_mesh = surface.mesh is not None
    coords = np.asarray(surface.mesh.vertices if is_mesh else surface.points_mm, dtype=float)
    faces = np.asarray(surface.mesh.faces, dtype=int) if is_mesh else None
    modes = [np.asarray(calculation.total_gy).ravel()] + [
        np.asarray(x).ravel() for x in calculation.event_maps_gy]
    maximum = float(np.max(modes[0]))
    if maximum <= 0 or not np.all(np.isfinite(coords)) or not all(np.all(np.isfinite(x)) for x in modes):
        raise ValueError("Visor 3D: coordenadas o dosis no finitas.")
    color_map = matplotlib.colormaps["inferno"]
    palette = ["#{:02x}{:02x}{:02x}".format(*np.rint(np.asarray(color_map(i/255)[:3])*255).astype(int))
               for i in range(256)]
    no_dose = "#505d70"

    def colors(dose: np.ndarray) -> list[str]:
        levels = np.clip((dose/maximum*255).astype(int), 0, 255)
        return [palette[int(level)] if value > 0 else no_dose
                for level, value in zip(levels, dose)]

    peaks = [surface.points_mm[int(np.argmax(m))].tolist() for m in modes]
    first_color = colors(modes[0])
    if is_mesh:
        fig = go.Figure(go.Mesh3d(
            x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
            i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
            facecolor=first_color, flatshading=True, hoverinfo="skip",
            lighting=dict(ambient=.82, diffuse=.65, specular=.12, roughness=.9),
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
    fig.update_layout(
        paper_bgcolor="#111722", font=dict(color="#e7eaf2", family="Arial, sans-serif"),
        margin=dict(l=0, r=0, t=0, b=0), showlegend=False,
        scene=dict(xaxis=dict(title="x lateral [mm]", backgroundcolor="#182131", color="#c6d0df"),
                   yaxis=dict(title="y craneal [mm]", backgroundcolor="#182131", color="#c6d0df"),
                   zaxis=dict(title="z anterior [mm]", backgroundcolor="#182131", color="#c6d0df"),
                   aspectmode="data", dragmode="orbit",
                   camera=dict(eye=dict(x=.35, y=.20, z=-2.25), up=dict(x=0, y=1, z=0))),
        uirevision="dose-map")
    figure = to_html(fig, full_html=False, include_plotlyjs=True, div_id="dose-map",
                     config=dict(displaylogo=False, responsive=True,
                                 modeBarButtonsToRemove=["lasso2d", "select2d"]))
    payload = {
        "kind": surface.kind,
        "vertices": np.round(coords, 3).tolist(),
        "faces": faces.tolist() if is_mesh else None,
        "doses": [np.round(m, 11).tolist() for m in modes],
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
    }
    template = files("skindose_verifier").joinpath("data/visor_3d_template.html").read_text(encoding="utf-8")
    # JSON contains only finite numerical arrays and fixed text labels; escaping
    # angle brackets also avoids terminating the application/json script element.
    encoded = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).replace("<", "\\u003c")
    path.write_text(template.replace("__PLOT_HTML__", figure).replace("__VIEWER_DATA__", encoded),
                    encoding="utf-8")
