import numpy as np
import trimesh
from trimesh.ray.ray_triangle import ray_triangle_id
from skindose_verifier.raytrace import first_hits


def test_segment_boxes_preserve_first_intersections():
    mesh = trimesh.creation.icosphere(subdivisions=3, radius=50.)
    targets = mesh.triangles_center[::29]
    for source in (np.array([150.,130.,120.]), np.array([0.,0.,150.]), np.array([0.,0.,0.])):
        vectors = targets-source
        distances = np.linalg.norm(vectors,axis=1)
        directions = vectors/distances[:,None]
        ids, rays, locations = ray_triangle_id(mesh.triangles, np.broadcast_to(source, targets.shape),
            directions, triangles_normal=mesh.face_normals, tree=mesh.triangles_tree, multiple_hits=False)
        expected = np.full(len(targets), -1)
        finite = np.linalg.norm(locations-source,axis=1) <= distances[rays]+1e-7
        expected[rays[finite]] = ids[finite]
        actual, stats = first_hits(mesh, source, targets)
        np.testing.assert_array_equal(actual,expected)
        assert stats['rays_traced'] == len(targets)


def test_rays_ending_before_object_do_not_hit_beyond_target():
    mesh = trimesh.creation.icosphere(subdivisions=2, radius=10.)
    actual, _ = first_hits(mesh, np.array([0.,0.,100.]), np.array([[0.,0.,50.],[100.,100.,50.]]))
    np.testing.assert_array_equal(actual, [-1,-1])
