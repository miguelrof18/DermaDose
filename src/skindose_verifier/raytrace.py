"""Instrumented Trimesh R-tree intersection, in bounded batches.

Counts are candidate ray/triangle pairs actually passed to the intersection
kernel, including candidates rejected by the triangle or forward-ray tests.
"""
import numpy as np
from trimesh.ray.ray_triangle import ray_triangle_id


class CountingTree:
    def __init__(self, tree):
        self.tree = tree
        self.bounds = tree.bounds
        self.candidate_pairs = 0

    def intersection(self, bounds):
        candidates = list(self.tree.intersection(bounds))
        self.candidate_pairs += len(candidates)
        return candidates


def first_hits(mesh, source_mm: np.ndarray, targets_mm: np.ndarray,
               batch_size: int = 64) -> tuple[np.ndarray, dict]:
    tree = CountingTree(mesh.triangles_tree)
    result = np.full(len(targets_mm), -1, dtype=int)
    for start in range(0, len(targets_mm), batch_size):
        direction = targets_mm[start:start+batch_size]-source_mm
        direction /= np.linalg.norm(direction, axis=1)[:, None]
        tri, ray, _ = ray_triangle_id(mesh.triangles, np.broadcast_to(source_mm, direction.shape),
                                     direction, triangles_normal=mesh.face_normals,
                                     tree=tree, multiple_hits=False)
        result[start+ray] = tri
    return result, {"rays_traced": len(targets_mm), "candidate_triangle_pairs": tree.candidate_pairs,
                    "triangles_in_mesh": len(mesh.faces), "backend": "trimesh_rtree_cpu"}
