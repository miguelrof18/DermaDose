"""Native VTK geometry checks; no Slicer GUI required."""
import sys
from pathlib import Path
import pytest
vtk = pytest.importorskip("vtk")
sys.path.insert(0,str(Path(__file__).resolve().parents[1] /
    "DermaDose_Surface/DermaDose_Surface/slicer_extension/DermaDoseSurface"))
from surface_processing import simplify_surface, surface_quality, write_stl_checked

def sphere(theta=450, phi=400):
    s=vtk.vtkSphereSource()
    s.SetRadius(100)
    s.SetThetaResolution(theta)
    s.SetPhiResolution(phi)
    s.Update()
    return s.GetOutput()

def test_range_geometry_and_disk(tmp_path):
    source=sphere()
    reduced, info=simplify_surface(source)
    assert info["export_allowed"], info
    assert 100000 <= info["output_triangles"] <= 200000
    assert info["max_vertex_surface_distance_mm"] <= 1
    assert info["relative_volume_change"] <= .01
    q=surface_quality(reduced)
    assert reduced.GetNumberOfPoints()-q["triangles"]//2 == 2
    disk=write_stl_checked(reduced,tmp_path/"surface.stl")
    assert disk["triangles"]==info["output_triangles"]

def test_small_not_subdivided():
    source=sphere(20,20)
    result, info=simplify_surface(source)
    assert info["below_range_without_subdivision"]
    assert result.GetNumberOfPolys()==source.GetNumberOfPolys()

def test_open_surface_rejected():
    plane=vtk.vtkPlaneSource()
    plane.Update()
    with pytest.raises(ValueError):
        simplify_surface(plane.GetOutput())

def test_geometric_failure_does_not_export_original():
    _,info=simplify_surface(sphere(),max_deviation_mm=1e-12)
    assert not info["export_allowed"]
    assert info["fallback_original"]
