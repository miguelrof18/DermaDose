import numpy as np
import pytest
import trimesh
from skindose_verifier.raytrace import first_hits
from skindose_verifier.vtk_ray_backend import VTKWorker, find_python_slicer


def test_optional_vtk_matches_trimesh_on_axial_oblique_and_interior_rays():
    executable = find_python_slicer()
    if executable is None:
        pytest.skip('PythonSlicer no instalado; se usa RTree')
    mesh = trimesh.creation.icosphere(subdivisions=3, radius=50.)
    worker = VTKWorker(mesh, executable)
    try:
        for source in (np.array([100.,120.,150.]),np.array([0.,0.,100.]),np.zeros(3)):
            targets = mesh.triangles_center[::17]
            expected, _ = first_hits(mesh, source, targets)
            np.testing.assert_array_equal(worker.intersect(source, targets), expected)
        actual = worker.intersect(np.array([0.,0.,100.]),np.array([[0.,0.,70.],[100.,100.,50.]]))
        np.testing.assert_array_equal(actual,[-1,-1])
    finally:
        worker.close()
