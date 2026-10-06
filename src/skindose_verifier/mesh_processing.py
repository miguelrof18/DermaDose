"""Conforming edge refinement without moving the supplied patient surface.

Splitting a shared edge on BOTH incident faces avoids T junctions. This is
linear subdivision, not an anatomical reconstruction or a smoothing filter.
See RESEARCH_REVIEW.md: Rana 2012; Trimesh remeshing documentation.
"""
from __future__ import annotations

import numpy as np
import trimesh

from .domain import InputError


def refine_mesh(mesh: trimesh.Trimesh, max_edge_mm: float, *,
                roi_mm: list[list[float]] | None = None,
                max_faces: int = 1_000_000) -> trimesh.Trimesh:
    """Bisect long shared edges, with optional axis-aligned refinement region.

Faces whose bounding boxes intersect the ROI drive refinement. Their neighbours
share the new vertices. The limit applies before allocation of the new faces.
"""
    if not np.isfinite(max_edge_mm) or max_edge_mm <= 0:
        raise InputError("mesh_max_edge_mm debe ser positivo y finito.")
    if not isinstance(max_faces, int) or max_faces < len(mesh.faces):
        raise InputError("mesh_max_faces debe admitir al menos la malla original.")
    roi = None if roi_mm is None else np.asarray(roi_mm, dtype=float)
    if roi is not None and (roi.shape != (2, 3) or not np.all(np.isfinite(roi)) or
                            np.any(roi[0] >= roi[1])):
        raise InputError("mesh_refinement_roi_mm requiere dos esquinas [mínimo, máximo].")
    vertices, faces = np.array(mesh.vertices), np.array(mesh.faces)
    parent = np.arange(len(faces))
    initial_count = len(faces)
    for iteration in range(24):
        edges, inverse = np.unique(np.sort(faces[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1),
                                   axis=0, return_inverse=True)
        inverse = inverse.reshape(-1, 3)
        lengths = np.linalg.norm(vertices[edges[:, 1]]-vertices[edges[:, 0]], axis=1)
        marked = lengths > max_edge_mm*(1+1e-12)
        if roi is not None:
            triangles = vertices[faces]
            eligible = np.all(triangles.max(axis=1) >= roi[0], axis=1) & np.all(
                triangles.min(axis=1) <= roi[1], axis=1)
            requested = np.zeros(len(edges), dtype=bool)
            requested[inverse[eligible].ravel()] = True
            marked &= requested
        if not marked.any():
            break
        flags = marked[inverse]
        counts = flags.sum(axis=1)
        next_count = int(np.sum(counts+1))
        if next_count > max_faces:
            raise InputError(f"Refinamiento excedería mesh_max_faces={max_faces}; aumentar paso o limitar región.")
        midpoint = np.full(len(edges), -1, dtype=int)
        midpoint[marked] = np.arange(len(vertices), len(vertices)+int(marked.sum()))
        vertices = np.vstack((vertices, vertices[edges[marked]].mean(axis=1)))
        mids = midpoint[inverse]
        blocks, ancestors = [faces[counts == 0]], [parent[counts == 0]]
        # Rotate each 1/2-edge pattern to a common winding-preserving template.
        for n in (1, 2):
            for start in range(3):
                mask = (counts == n) & flags[:, start]
                mask &= (~flags[:, (start+1) % 3] if n == 1 else flags[:, (start+1) % 3])
                f, m = np.roll(faces[mask], -start, axis=1), np.roll(mids[mask], -start, axis=1)
                a, b, c = f.T
                ab, bc, _ = m.T
                if n == 1:
                    block = np.stack((np.column_stack((a, ab, c)), np.column_stack((ab, b, c))), axis=1)
                else:
                    block = np.stack((np.column_stack((a, ab, c)), np.column_stack((ab, bc, c)),
                                      np.column_stack((ab, b, bc))), axis=1)
                blocks.append(block.reshape(-1, 3)); ancestors.append(np.repeat(parent[mask], n+1))
        mask = counts == 3
        a, b, c = faces[mask].T
        ab, bc, ca = mids[mask].T
        blocks.append(np.stack((np.column_stack((a, ab, ca)), np.column_stack((ab, b, bc)),
                                np.column_stack((ca, bc, c)), np.column_stack((ab, bc, ca))), axis=1).reshape(-1, 3))
        ancestors.append(np.repeat(parent[mask], 4))
        faces, parent = np.vstack(blocks), np.concatenate(ancestors)
    else:
        raise InputError("Refinamiento no convergió en 24 iteraciones.")
    refined = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    if not refined.is_watertight or not refined.is_winding_consistent or not refined.is_volume:
        raise InputError("El refinamiento no conservó una superficie cerrada y orientada.")
    if not np.isclose(refined.area, mesh.area, rtol=1e-10) or not np.isclose(refined.volume, mesh.volume, rtol=1e-10):
        raise InputError("El refinamiento alteró el área o el volumen del fantoma.")
    refined.metadata.update(original_face_count=initial_count, refinement_iterations=iteration,
                            max_edge_target_mm=max_edge_mm, refinement_roi_mm=roi_mm,
                            parent_face_index=parent)
    return refined


def mesh_quality(mesh: trimesh.Trimesh) -> dict:
    """Diagnostics, not a certificate that the mesh is free of intersections."""
    edges = mesh.edges_unique_length
    angles = np.degrees(mesh.face_angles)
    return {"face_count": int(len(mesh.faces)), "vertex_count": int(len(mesh.vertices)),
            "watertight": bool(mesh.is_watertight), "winding_consistent": bool(mesh.is_winding_consistent),
            "area_cm2": float(mesh.area/100), "volume_litre": float(mesh.volume/1e6),
            "edge_p95_mm": float(np.quantile(edges, .95)), "maximum_edge_mm": float(edges.max()),
            "minimum_triangle_angle_deg": float(angles.min()),
            "triangle_area_p95_mm2": float(np.quantile(mesh.area_faces, .95)),
            "maximum_triangle_area_mm2": float(mesh.area_faces.max()),
            "bounds_patient_mm": mesh.bounds.tolist(),
            "original_face_count": int(mesh.metadata.get("original_face_count", len(mesh.faces))),
            "max_edge_target_mm": mesh.metadata.get("max_edge_target_mm"),
            "refinement_roi_mm": mesh.metadata.get("refinement_roi_mm"),
            "self_intersections_checked": False,
            "anatomical_smoothing_applied": False}
