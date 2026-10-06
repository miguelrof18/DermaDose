"""Informes anonimizados y mapas. Nunca exporta PatientName/PatientID."""
from __future__ import annotations

from pathlib import Path
import csv
import json

import numpy as np

from .engine import Calculation


def save_results(calculation: Calculation, directory: str | Path, *, audit: dict) -> dict[str, str]:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    report = {**calculation.summary, "audit": audit}
    json_path = directory/"resultados.json"
    json_path.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    raw_path = directory/"mapa_dosis.npz"
    np.savez_compressed(raw_path, dose_gy=calculation.total_gy.reshape(calculation.surface.shape),
                        event_dose_gy=np.stack(calculation.event_maps_gy).reshape(
                            (len(calculation.event_maps_gy),) + calculation.surface.shape),
                        event_uid=np.asarray([ev["event_uid"] for ev in calculation.summary["events"]]),
                        patient_points_mm=calculation.surface.points_mm,
                        normals=calculation.surface.normals, area_mm2=calculation.surface.area_mm2,
                        axis_rows=calculation.surface.axes[0], axis_columns=calculation.surface.axes[1],
                        vertices_mm=(calculation.surface.mesh.vertices if calculation.surface.mesh is not None else np.zeros((0, 3))),
                        faces=(calculation.surface.mesh.faces if calculation.surface.mesh is not None else np.zeros((0, 3), dtype=int)))
    csv_path = directory/"eventos.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["event_index", "plane", "ka_r_gy", "kap_gy_m2",
                                                    "maximum_gy", "irradiated_cells", "sweep_assumed"])
        writer.writeheader()
        for idx, ev in enumerate(calculation.summary["events"], start=1):
            writer.writerow({"event_index": idx, "plane": ev["plane"], "ka_r_gy": ev["ka_r_gy"],
                             "kap_gy_m2": ev["kap_gy_m2"], "maximum_gy": ev["maximum_gy"],
                             "irradiated_cells": ev["irradiated_cells"],
                             "sweep_assumed": ev["uniform_sweep_assumed"]})
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, ListedColormap
    from matplotlib.cm import ScalarMappable
    if calculation.surface.kind == "mesh":
        shades = plt.colormaps["inferno"](np.linspace(0, 1, 256))
        shades[0] = [.26, .29, .32, 1.]
        palette = ListedColormap(shades)
    else:
        palette = plt.colormaps["inferno"]
    fig, ax = plt.subplots(figsize=(6, 8) if calculation.surface.kind == "mesh" else (10, 5))
    ax.set_title(f"{calculation.summary['quantity']} (INVESTIGACIÓN, SIN VALIDACIÓN CLÍNICA)")
    dose = calculation.total_gy.reshape(calculation.surface.shape)
    if calculation.surface.kind == "mesh":
        points = calculation.surface.points_mm
        # A projection overlays surfaces; the 3D picture contains the real faces.
        order = np.argsort(calculation.total_gy)
        im = ax.scatter(points[order, 0], points[order, 1], c=calculation.total_gy[order],
                        s=3, cmap=palette, vmin=0, linewidths=0)
        ax.set_aspect("equal")
        ax.set_xlabel("x lateral paciente [mm]")
        ax.set_title("Proyección de caras triangulares\n(superficies superpuestas)", fontsize=11)
    else:
        y, x = calculation.surface.axes
        extent = [float(x.min()), float(x.max()), float(y.min()), float(y.max())]
        im = ax.imshow(dose, origin="lower", aspect="auto", extent=extent, cmap=palette,
                       interpolation="nearest", vmin=0)
        if calculation.surface.kind == "ellipse":
            ax.set_xlabel("ángulo circunferencial θ [rad]")
        else:
            ax.set_xlabel("posición lateral local [mm]")
    ax.set_ylabel("posición longitudinal local [mm]")
    fig.colorbar(im, ax=ax, label="Gy")
    fig.tight_layout()
    map_path = directory/"mapa_2d.png"
    fig.savefig(map_path, dpi=150)
    plt.close(fig)
    # Three-dimensional spatial depiction uses the same dose values/point order.
    fig = plt.figure(figsize=(8, 7))
    ax3 = fig.add_subplot(111, projection="3d")
    points = calculation.surface.points_mm
    norm = Normalize(vmin=0, vmax=max(1e-12, float(np.max(calculation.total_gy))))
    if calculation.surface.mesh is not None:
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        mesh = calculation.surface.mesh
        collection = Poly3DCollection(np.asarray(mesh.vertices)[mesh.faces],
                                      facecolors=palette(norm(calculation.total_gy)),
                                      edgecolors="none", linewidths=0)
        ax3.add_collection3d(collection)
        ax3.auto_scale_xyz(*mesh.bounds.T)
        ax3.set_box_aspect(np.ptp(np.asarray(mesh.vertices), axis=0))
        ax3.view_init(elev=-65, azim=90)
        sc = ScalarMappable(norm=norm, cmap=palette)
    else:
        stride = max(1, len(points)//15000)
        sc = ax3.scatter(points[::stride, 0], points[::stride, 1], points[::stride, 2],
                         c=calculation.total_gy[::stride], s=4, cmap=palette, norm=norm)
    ax3.set_xlabel("x [mm]")
    ax3.set_ylabel("y [mm]")
    ax3.set_zlabel("z [mm]")
    fig.colorbar(sc, ax=ax3, label="Gy", shrink=.7)
    fig.suptitle("Distribución de dosis, geometría de investigación")
    fig.tight_layout()
    map3_path = directory/"mapa_3d.png"
    fig.savefig(map3_path, dpi=150)
    plt.close(fig)
    from .interactive import save_interactive
    histogram = calculation.summary["surface_metrics"]["cumulative_dose_area_histogram"]
    dah_path = directory/"histograma_dosis_area.csv"
    with dah_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["dose_threshold_gy", "area_at_least_cm2"])
        writer.writerows(zip(histogram["dose_gy"], histogram["area_at_least_cm2"]))
    fig, axis = plt.subplots(figsize=(7, 4))
    areas = list(histogram["area_at_least_cm2"])
    # At 0+, restrict plotting to irradiated area. CSV retains A(D>=0).
    areas[0] = calculation.summary["surface_metrics"]["irradiated_surface"]["area_cm2"]
    axis.step(np.asarray(histogram["dose_gy"])*1000, areas, where="post", color="#b95a09")
    axis.set(xlabel="Umbral [mGy]", ylabel="Superficie [cm²]", title="Histograma dosis–área · piel irradiada")
    axis.grid(alpha=.2); fig.tight_layout()
    dah_image = directory/"histograma_dosis_area.png"
    fig.savefig(dah_image, dpi=150); plt.close(fig)
    interactive_path = directory/"mapa_3d_interactivo.html"
    save_interactive(calculation, interactive_path)
    return {k: str(v) for k, v in {"summary": json_path, "numeric_map": raw_path,
                                     "events": csv_path, "map_2d": map_path, "map_3d": map3_path,
                                     "interactive_3d": interactive_path, "dose_area_csv": dah_path,
                                     "dose_area_plot": dah_image}.items()}
