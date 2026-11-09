"""Instrumented Trimesh R-tree intersection, in bounded batches.

Counts are candidate ray/triangle pairs actually passed to the intersection
kernel, including candidates rejected by the triangle or forward-ray tests.
"""
import numpy as np
from trimesh.ray.ray_triangle import ray_triangle_id


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
