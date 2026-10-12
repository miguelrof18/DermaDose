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

        candidates = self.tree.intersection(tuple(np.r_[query_min, query_max]))

        def counted_candidates():
            for candidate in candidates:
                self.candidate_pairs += 1
                yield candidate

        return counted_candidates()


def first_hits(mesh, source_mm: np.ndarray, targets_mm: np.ndarray,
               batch_size: int = 64) -> tuple[np.ndarray, dict]:
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
