"""Static orthographic views of a triangulated skin surface with isodose curves.

Each view is rasterised with a depth buffer, so the faces nearest to the
observer hide the ones behind them (a posterior field is no longer covered by
the anterior skin, as happens with a 3D point cloud). Only the picture is
built here: the dose per face is the computed value, never resampled into the
results.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

ISODOSE_PERCENT = (25.0, 50.0, 75.0)
# Dark/white pairs keep each curve readable over green, yellow, orange and red.
ISODOSE_STYLE = (("#1d3fb8", "#ffffff"), ("#0b0f17", "#ffffff"), ("#ffffff", "#0b0f17"))


MIN_ISLAND_CM2 = 2.0


def _clean_region(region: np.ndarray, valid: np.ndarray, pixel_area_cm2: float) -> np.ndarray:
    """Region at or above an isodose level, for drawing only.

    Pinholes left by mesh folds and islands smaller than MIN_ISLAND_CM2 are
    removed so each curve follows the field edge instead of every fold.
    """
    from scipy import ndimage
    min_pixels = max(1, int(round(MIN_ISLAND_CM2/pixel_area_cm2)))
    region = ndimage.binary_closing(region, structure=np.ones((3, 3)), iterations=2) & valid
    labels, count = ndimage.label(region)
    if count:
        sizes = ndimage.sum(region, labels, index=np.arange(1, count + 1))
        region = np.isin(labels, 1 + np.flatnonzero(sizes >= min_pixels))
    holes = ~region & valid
    labels, count = ndimage.label(holes)
    if count:
        sizes = ndimage.sum(holes, labels, index=np.arange(1, count + 1))
        region |= np.isin(labels, 1 + np.flatnonzero(sizes < min_pixels))
    return region


def _smooth_mask(region: np.ndarray) -> np.ndarray:
    from scipy import ndimage
    return ndimage.gaussian_filter(region.astype(float), sigma=1.2)


@dataclass
class View:
    name: str
    image_face: np.ndarray       # face index per pixel, -1 outside the surface
    extent: tuple[float, float, float, float]
    left_label: str
    right_label: str
    forward: np.ndarray


def _anatomy(ct_lps: bool) -> tuple[np.ndarray, np.ndarray, bool]:
    """Superior and posterior unit vectors, and whether left/right are known."""
    if ct_lps:
        return np.array([0., 0., 1.]), np.array([0., 1., 0.]), True
    # Engine axes of the provisional surface: head towards +y, anterior +z.
    return np.array([0., 1., 0.]), np.array([0., 0., -1.]), False


def _rasterise(vertices: np.ndarray, faces: np.ndarray, forward: np.ndarray, up: np.ndarray,
               pixel_mm: float) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    right = np.cross(forward, up)
    u, v, depth = vertices @ right, vertices @ up, vertices @ forward
    u0, v0 = float(u.min()) - pixel_mm, float(v.min()) - pixel_mm
    width = int(np.ceil((u.max() - u0)/pixel_mm)) + 2
    height = int(np.ceil((v.max() - v0)/pixel_mm)) + 2
    uv = np.stack([u, v], axis=1)
    tri = uv[faces]
    edge = np.max(np.linalg.norm(tri - np.roll(tri, 1, axis=1), axis=2), axis=1)
    # Barycentric sub-samples per face, dense enough to leave no gaps.
    subdivisions = np.clip(np.ceil(edge/pixel_mm*1.5).astype(int), 1, 12)
    pixel_ids, depths, face_ids = [], [], []
    for n in np.unique(subdivisions):
        selected = np.flatnonzero(subdivisions == n)
        grid = np.array([(i, j) for i in range(n + 1) for j in range(n + 1 - i)], dtype=float)/n
        weights = np.column_stack([grid, 1 - grid.sum(axis=1)])            # (s, 3)
        points_u = weights @ u[faces[selected]].T                          # (s, f)
        points_v = weights @ v[faces[selected]].T
        points_d = weights @ depth[faces[selected]].T
        cols = np.clip(((points_u - u0)/pixel_mm).astype(int), 0, width - 1)
        rows = np.clip(((points_v - v0)/pixel_mm).astype(int), 0, height - 1)
        pixel_ids.append((rows*width + cols).ravel())
        depths.append(points_d.ravel())
        face_ids.append(np.broadcast_to(selected, points_u.shape).ravel())
    pixel_ids, depths, face_ids = map(np.concatenate, (pixel_ids, depths, face_ids))
    order = np.lexsort((depths, pixel_ids))       # nearest sample first per pixel
    first = np.ones(order.size, dtype=bool)
    first[1:] = pixel_ids[order][1:] != pixel_ids[order][:-1]
    image = np.full(width*height, -1, dtype=np.int64)
    image[pixel_ids[order][first]] = face_ids[order][first]
    image = image.reshape(height, width)
    # Close pinholes left between samples, only inside the silhouette.
    from scipy import ndimage
    covered = image >= 0
    inside = ndimage.binary_closing(covered, structure=np.ones((3, 3)), iterations=2) | covered
    holes = inside & ~covered
    if holes.any():
        _, (ri, ci) = ndimage.distance_transform_edt(~covered, return_indices=True)
        image[holes] = image[ri[holes], ci[holes]]
    extent = (u0, u0 + width*pixel_mm, v0, v0 + height*pixel_mm)
    return image, extent


def build_views(vertices: np.ndarray, faces: np.ndarray, *, ct_lps: bool,
                pixel_mm: float = 1.5) -> list[View]:
    vertices = np.asarray(vertices, dtype=float)
    faces = np.asarray(faces, dtype=int)
    superior, posterior, sided = _anatomy(ct_lps)
    left = np.cross(posterior, superior)
    lateral = ("Lateral izquierda", "Lateral derecha") if sided else ("Lateral +x", "Lateral −x")
    views = []
    for name, forward in (("Posterior", -posterior), ("Anterior", posterior),
                          (lateral[0], -left), (lateral[1], left)):
        image, extent = _rasterise(vertices, faces, forward, superior, pixel_mm)
        right = np.cross(forward, superior)
        if abs(right @ left) > .5:
            labels = ("D", "I") if right @ left > 0 else ("I", "D")
            if not sided:
                labels = ("−x", "+x") if right @ left > 0 else ("+x", "−x")
        else:
            labels = ("A", "P") if right @ posterior > 0 else ("P", "A")
        views.append(View(name, image, extent, *labels, forward))
    return views


def save_surface_views(path, vertices, faces, face_dose_gy, face_normals, *, ct_lps: bool,
                       palette, no_dose_rgb=(.26, .29, .32), title: str = "") -> None:
    import matplotlib.pyplot as plt
    from matplotlib import patheffects
    from matplotlib.lines import Line2D
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    dose = np.asarray(face_dose_gy, dtype=float).ravel()
    normals = np.asarray(face_normals, dtype=float)
    peak = float(dose.max())
    norm = Normalize(vmin=0, vmax=max(peak, 1e-12))
    views = build_views(vertices, faces, ct_lps=ct_lps)
    widths = [v.extent[1] - v.extent[0] for v in views]
    fig, axes = plt.subplots(1, len(views), figsize=(15, 7.2),
                             gridspec_kw=dict(width_ratios=widths))
    levels = [peak*p/100 for p in ISODOSE_PERCENT] if peak > 0 else []
    for index, (ax, view) in enumerate(zip(axes, views)):
        forward = view.forward
        face = view.image_face
        valid = face >= 0
        safe = np.where(valid, face, 0)
        pixel_dose = np.where(valid, dose[safe], np.nan)
        rgb = np.asarray(palette(norm(np.nan_to_num(pixel_dose))))[..., :3]
        rgb[valid & (pixel_dose <= 0)] = no_dose_rgb
        # Light shading by obliquity gives the body its shape without hiding colour.
        shade = .62 + .38*np.abs(normals[safe] @ forward)
        rgba = np.ones(face.shape + (4,))
        rgba[..., :3] = np.clip(rgb*shade[..., None], 0, 1)
        rgba[~valid] = (1, 1, 1, 0)
        ax.imshow(rgba, origin="lower", extent=view.extent, interpolation="nearest")
        peak_face = int(np.argmax(dose))
        peak_xy = None
        if peak > 0 and np.any(face == peak_face):
            r, c = np.argwhere(face == peak_face)[0]
            peak_xy = np.array([view.extent[0] + (c + .5)*(view.extent[1]-view.extent[0])/face.shape[1],
                                view.extent[2] + (r + .5)*(view.extent[3]-view.extent[2])/face.shape[0]])
        if levels and np.nanmax(pixel_dose) >= levels[0]:
            height, width = face.shape
            xs = np.linspace(view.extent[0], view.extent[1], width, endpoint=False) + (view.extent[1]-view.extent[0])/width/2
            ys = np.linspace(view.extent[2], view.extent[3], height, endpoint=False) + (view.extent[3]-view.extent[2])/height/2
            dose_px = np.where(valid, np.nan_to_num(pixel_dose), 0.)
            pixel_area_cm2 = ((view.extent[1]-view.extent[0])/width)*((view.extent[3]-view.extent[2])/height)/100
            for level, percent, (line, halo) in zip(levels, ISODOSE_PERCENT, ISODOSE_STYLE):
                region = _clean_region(dose_px >= level, valid, pixel_area_cm2)
                if not region.any():
                    continue
                contour = ax.contour(xs, ys, _smooth_mask(region), levels=[.5], colors=[line],
                                     linewidths=1.3)
                for artist in getattr(contour, "collections", [contour]):
                    artist.set_path_effects([patheffects.Stroke(linewidth=2.6, foreground=halo, alpha=.85),
                                             patheffects.Normal()])
                columns = np.flatnonzero(region.any(axis=0))
                if columns.size*(view.extent[1]-view.extent[0])/width < 60:
                    continue  # too narrow for a readable label (seen edge-on)
                labels = ax.clabel(contour, fmt={.5: f"{percent:g} %"}, fontsize=7.5, inline=True,
                                   inline_spacing=2, colors=[line], manual=False)
                for label in list(labels):
                    if peak_xy is not None and np.hypot(*(np.asarray(label.get_position()) - peak_xy)) < 18:
                        label.remove()
                        continue
                    label.set_path_effects([patheffects.withStroke(linewidth=2.2, foreground=halo)])
        if peak_xy is not None:
            ax.plot([peak_xy[0]], [peak_xy[1]], marker="o", ms=7, mfc="white", mec="#0b0f17", mew=1.6)
        ax.set_title(view.name, fontsize=11)
        ax.set_xticks([])
        ax.set_aspect("equal")
        for side in ("top", "right", "left", "bottom"):
            ax.spines[side].set_visible(False)
        ax.text(0, -.02, view.left_label, transform=ax.transAxes, ha="left", va="top", fontsize=10, color="#374151")
        ax.text(1, -.02, view.right_label, transform=ax.transAxes, ha="right", va="top", fontsize=10, color="#374151")
        if index == 0:
            ax.set_ylabel("z craneal TC [mm]" if ct_lps else "y craneal [mm]")
        else:
            ax.set_yticks([])
    handles = [Line2D([], [], color=line, lw=1.3,
                      path_effects=[patheffects.Stroke(linewidth=2.6, foreground=halo), patheffects.Normal()],
                      label=f"{p:g} % ({level*1000:.0f} mGy)")
               for p, level, (line, halo) in zip(ISODOSE_PERCENT, levels, ISODOSE_STYLE)]
    handles.append(Line2D([], [], ls="", marker="o", ms=7, mfc="white", mec="#0b0f17", mew=1.6,
                          label=f"Máximo ({peak*1000:.0f} mGy)"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False,
               title="Isodosis respecto al máximo de la dosis acumulada", fontsize=9, title_fontsize=9)
    fig.subplots_adjust(left=.05, right=.9, bottom=.14, top=.9, wspace=.06)
    sm = ScalarMappable(norm=norm, cmap=palette)
    fig.colorbar(sm, cax=fig.add_axes((.925, .22, .012, .6)), label="Gy")
    if title:
        fig.suptitle(title)
    fig.savefig(path, dpi=150)
    plt.close(fig)
