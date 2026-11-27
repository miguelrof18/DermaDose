"""First-entry visibility on closed meshes.

Two kernels are available:

* ``visible_faces``: vectorised central projection from the focal spot. For a
  closed, consistently oriented surface and a source outside it, the first
  crossing of a source-to-target segment is always an *entry* through a face
  that looks at the source. A target face is therefore visible exactly when no
  other source-facing triangle is crossed before it. The projection through the
  focal spot preserves ray/triangle incidence, so the test reduces to a 2D
  point-in-triangle query plus a depth comparison, evaluated in NumPy on a
  uniform grid of the projected field. No physics is added or removed.
* ``first_hits``: the instrumented Trimesh R-tree path, kept for arbitrary
  targets, as reference implementation and as automatic fallback.

Counts are candidate ray/triangle pairs actually passed to the intersection
kernel, including candidates rejected by the triangle or forward-ray tests.
"""
import os

import numpy as np


def _beam_basis(axis):
    axis = axis/np.linalg.norm(axis)
    helper = np.array([1., 0., 0.]) if abs(axis[0]) < .9 else np.array([0., 1., 0.])
    first = np.cross(axis, helper)
    first /= np.linalg.norm(first)
    return first, np.cross(axis, first), axis


def visible_faces(vertices, faces, face_normals, source_mm, target_faces, *,
                  max_pairs_per_chunk=4_000_000):
    """Return (visible_mask, statistics) for faces ``target_faces``, or None.

    ``None`` means that the projective test is not applicable (a source-facing
    triangle reaches the plane of the focal spot); callers then use the
    reference ray tracer. Targets are the face centroids.
    """
    vertices = np.asarray(vertices, dtype=float)
    faces = np.asarray(faces, dtype=np.int64)
    normals = np.asarray(face_normals, dtype=float)
    targets = np.asarray(target_faces, dtype=np.int64)
    source = np.asarray(source_mm, dtype=float)
    stats = {"rays_traced": int(len(targets)), "candidate_triangle_pairs": 0,
             "candidate_triangle_pairs_available": True, "triangles_in_mesh": int(len(faces)),
             "backend": "numpy_central_projection_cpu"}
    if not len(targets):
        return np.zeros(0, dtype=bool), stats
    relative = vertices - source
    corner0 = relative[faces[:, 0]]
    facing = np.einsum("ij,ij->i", normals, corner0) < 0.
    centroid = relative[faces[targets]].mean(axis=1)
    axis = centroid.mean(axis=0)
    if not np.all(np.isfinite(axis)) or np.linalg.norm(axis) == 0:
        return None
    a, b, w = _beam_basis(axis)
    depth_v = relative @ w
    target_depth = centroid @ w
    if np.any(target_depth <= 0):
        return None
    tri = np.flatnonzero(facing)
    tri_depth = depth_v[faces[tri]]
    scale = float(np.max(np.abs(depth_v)))
    eps_depth = 1e-9*scale
    in_front = tri_depth.max(axis=1) > eps_depth
    if np.any(in_front & (tri_depth.min(axis=1) <= eps_depth)):
        return None  # a source-facing triangle straddles the focal plane
    tri = tri[in_front & (tri_depth.min(axis=1) < target_depth.max())]
    # Central projection of vertices and targets onto the plane w = 1.
    tu = (centroid @ a)/target_depth
    tv = (centroid @ b)/target_depth
    visible = np.ones(len(targets), dtype=bool)
    if not len(tri):
        return visible, stats
    corners = faces[tri]
    pu = (relative[corners] @ a)/depth_v[corners]
    pv = (relative[corners] @ b)/depth_v[corners]
    lo_u, hi_u, lo_v, hi_v = pu.min(1), pu.max(1), pv.min(1), pv.max(1)
    span_u = max(float(tu.max()-tu.min()), 1e-12)
    span_v = max(float(tv.max()-tv.min()), 1e-12)
    pad = 1e-9*max(span_u, span_v, 1e-6)
    keep = ((hi_u >= tu.min()-pad) & (lo_u <= tu.max()+pad) &
            (hi_v >= tv.min()-pad) & (lo_v <= tv.max()+pad))
    tri, pu, pv = tri[keep], pu[keep], pv[keep]
    lo_u, hi_u, lo_v, hi_v = lo_u[keep], hi_u[keep], lo_v[keep], hi_v[keep]
    if not len(tri):
        return visible, stats
    # Uniform grid over the projected targets; cell ~ typical triangle size.
    extent = np.maximum(hi_u-lo_u, hi_v-lo_v)
    # Never more grid cells than about twice the number of targets: a large
    # triangle then expands into a bounded number of (cell, triangle) pairs.
    cell = max(float(np.median(extent))*1.5, np.sqrt(span_u*span_v/(2*len(targets)+16)), 1e-12)
    nu = int(min(max(1, np.ceil(span_u/cell)), 4096))
    nv = int(min(max(1, np.ceil(span_v/cell)), 4096))
    u0, v0 = float(tu.min()), float(tv.min())
    cu, cv = span_u/nu, span_v/nv

    def cells(values, origin, size, count):
        return np.clip(np.floor((values-origin)/size).astype(np.int64), 0, count-1)

    i0, i1 = cells(lo_u-pad, u0, cu, nu), cells(hi_u+pad, u0, cu, nu)
    j0, j1 = cells(lo_v-pad, v0, cv, nv), cells(hi_v+pad, v0, cv, nv)
    target_cell = cells(tu, u0, cu, nu)*nv + cells(tv, v0, cv, nv)
    # Expand each triangle into the cells covered by its projected box.
    width_i, width_j = i1-i0+1, j1-j0+1
    counts = width_i*width_j
    owner = np.repeat(np.arange(len(tri)), counts)
    local = np.arange(int(counts.sum())) - np.repeat(np.cumsum(counts)-counts, counts)
    cell_of_pair = (i0[owner] + local // width_j[owner])*nv + (j0[owner] + local % width_j[owner])
    order = np.argsort(cell_of_pair, kind="stable")
    cell_sorted, owner_sorted = cell_of_pair[order], owner[order]
    start = np.searchsorted(cell_sorted, target_cell, side="left")
    stop = np.searchsorted(cell_sorted, target_cell, side="right")
    lengths = stop-start
    # Normal and plane offset of each candidate, relative to the source.
    n = normals[tri]
    plane = np.einsum("ij,ij->i", n, relative[faces[tri, 0]])
    target_index = np.repeat(np.arange(len(targets)), lengths)
    offsets = np.arange(int(lengths.sum())) - np.repeat(np.cumsum(lengths)-lengths, lengths)
    candidate = owner_sorted[np.repeat(start, lengths) + offsets]
    stats["candidate_triangle_pairs"] = int(len(candidate))
    for begin in range(0, len(candidate), max_pairs_per_chunk):
        t_idx = target_index[begin:begin+max_pairs_per_chunk]
        c_idx = candidate[begin:begin+max_pairs_per_chunk]
        other = tri[c_idx] != targets[t_idx]
        t_idx, c_idx = t_idx[other], c_idx[other]
        x, y = tu[t_idx], tv[t_idx]
        ax, ay = pu[c_idx, 0], pv[c_idx, 0]
        bx, by = pu[c_idx, 1], pv[c_idx, 1]
        cx, cy = pu[c_idx, 2], pv[c_idx, 2]
        denominator = (by-cy)*(ax-cx) + (cx-bx)*(ay-cy)
        valid = np.abs(denominator) > 0
        safe = np.where(valid, denominator, 1.)
        l1 = ((by-cy)*(x-cx) + (cx-bx)*(y-cy))/safe
        l2 = ((cy-ay)*(x-cx) + (ax-cx)*(y-cy))/safe
        l3 = 1.-l1-l2
        tolerance = 1e-12
        inside = valid & (l1 >= -tolerance) & (l2 >= -tolerance) & (l3 >= -tolerance)
        # Segment parameter of the crossing with the triangle plane: < 1 means
        # the triangle lies between the focal spot and the target centroid.
        along = np.einsum("ij,ij->i", n[c_idx], centroid[t_idx])
        with np.errstate(divide="ignore", invalid="ignore"):
            parameter = plane[c_idx]/along
        hit = inside & (along < 0) & (parameter > 0) & (parameter < 1.-1e-7)
        visible[t_idx[hit]] = False
    return visible, stats


def _projective_enabled():
    return os.environ.get("DERMADOSE_RAY_BACKEND", "projective").strip().lower() == "projective"


class CountingTree:
    """Count candidates while restricting each query to its finite ray segment."""

    def __init__(self, tree, origins: np.ndarray, directions: np.ndarray,
                 max_distance: np.ndarray):
        self.tree = tree
        self.bounds = tree.bounds
        self.candidate_pairs = 0
        self.origins = origins
        self.directions = directions
        self.max_distance = max_distance
        self.query_index = 0

    def intersection(self, bounds):
        # Trimesh queries the tree once per ray. Its default ray bounds extend
        # through the entire mesh AABB; clipping them to the evaluated segment
        # avoids considering triangles beyond the surface sample.
        index = self.query_index
        self.query_index += 1
        origin = self.origins[index]
        end = origin + self.directions[index] * self.max_distance[index]
        segment_min = np.minimum(origin, end)
        segment_max = np.maximum(origin, end)
        query_min = np.maximum(np.asarray(bounds[:3]), segment_min)
        query_max = np.minimum(np.asarray(bounds[3:]), segment_max)
        if np.any(query_min > query_max):
            return ()

        # An oblique ray has a very large single AABB on dense CT meshes.
        # Cover the exact finite segment with smaller boxes; their union still
        # contains every possible hit, without changing the triangle kernel.
        direction = self.directions[index]
        mesh_bounds = np.asarray(self.bounds, dtype=float).reshape(2, 3)
        moving = np.abs(direction) > 1e-15
        if np.any(~moving & ((origin < mesh_bounds[0]) | (origin > mesh_bounds[1]))):
            return ()
        enter, leave = 0., float(self.max_distance[index])
        if np.any(moving):
            t = (mesh_bounds[:, moving] - origin[moving]) / direction[moving]
            enter = max(enter, float(np.min(t, axis=0).max()))
            leave = min(leave, float(np.max(t, axis=0).min()))
        if enter > leave:
            return ()
        segments = min(64, max(1, int(np.ceil((leave-enter)/10.))))
        points = origin + np.linspace(enter, leave, segments + 1)[:, None] * direction
        lows = np.maximum(np.minimum(points[:-1], points[1:]) - 1e-7, query_min)
        highs = np.minimum(np.maximum(points[:-1], points[1:]) + 1e-7, query_max)
        keep = np.all(lows <= highs, axis=1)
        lows, highs = lows[keep], highs[keep]
        if not len(lows):
            return ()
        if hasattr(self.tree, 'intersection_v'):
            try:
                ids, _ = self.tree.intersection_v(lows, highs)
                candidates = np.unique(ids)
            except (TypeError, AttributeError):
                candidates = np.unique([i for lo, hi in zip(lows, highs)
                                        for i in self.tree.intersection(tuple(np.r_[lo, hi]))])
        else:
            candidates = np.unique([i for lo, hi in zip(lows, highs)
                                    for i in self.tree.intersection(tuple(np.r_[lo, hi]))])
        self.candidate_pairs += len(candidates)
        return candidates


def first_hits(mesh, source_mm: np.ndarray, targets_mm: np.ndarray,
               batch_size: int = 64) -> tuple[np.ndarray, dict]:
    from trimesh.ray.ray_triangle import ray_triangle_id
    if len(mesh.faces) >= 200_000:
        from .vtk_ray_backend import intersect, close_worker
        try:
            result = intersect(mesh, source_mm, targets_mm)
            if result is not None:
                hits, version = result
                return hits, {"rays_traced": len(targets_mm), "candidate_triangle_pairs": 0,
                              "candidate_triangle_pairs_available": False,
                              "triangles_in_mesh": len(mesh.faces), "backend": "vtk_static_cell_locator",
                              "backend_version": version}
        except (OSError, ValueError, RuntimeError):
            close_worker()
    base_tree = mesh.triangles_tree
    result = np.full(len(targets_mm), -1, dtype=int)
    candidate_pairs = 0
    for start in range(0, len(targets_mm), batch_size):
        targets = targets_mm[start:start+batch_size]
        origins = np.broadcast_to(source_mm, targets.shape)
        vectors = targets-source_mm
        distances = np.linalg.norm(vectors, axis=1)
        direction = vectors / distances[:, None]
        tree = CountingTree(base_tree, origins, direction, distances)
        tri, ray, _ = ray_triangle_id(mesh.triangles, np.broadcast_to(source_mm, direction.shape),
                                     direction, triangles_normal=mesh.face_normals,
                                     tree=tree, multiple_hits=False)
        result[start+ray] = tri
        candidate_pairs += tree.candidate_pairs
    return result, {"rays_traced": len(targets_mm), "candidate_triangle_pairs": candidate_pairs,
                    "triangles_in_mesh": len(mesh.faces), "backend": "trimesh_rtree_cpu"}
