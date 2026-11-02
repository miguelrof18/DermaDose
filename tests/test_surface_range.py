"""Native VTK geometry checks; no Slicer GUI required."""
import sys
from pathlib import Path
import pytest
vtk = pytest.importorskip("vtk")
sys.path.insert(0,str(Path(__file__).resolve().parents[1] / "DermaDose_Surface/DermaDose_Surface/slicer_extension/DermaDoseSurface"))
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

def test_comparison_same_original_no_overwrite(tmp_path):
    from surface_pair import export_pair
    source=sphere()
    original=tmp_path/"original.vtp"
    writer=vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(original))
    writer.SetInputData(source)
    assert writer.Write()
    result=export_pair(original,tmp_path)
    low=result["meshes"]["surface_150k"]
    high=result["meshes"]["surface_300k"]
    assert 100000 <= low["output_triangles"] <= 200000
    assert 270000 <= high["output_triangles"] <= 330000
    assert low["source_sha256"] == high["source_sha256"]
    assert not result["psd_compared"]
    with pytest.raises(FileExistsError):
        export_pair(original,tmp_path)


def test_optimized_shape_and_candidate_selection():
    from surface_processing import triangle_statistics
    source = sphere(150,100)
    reduced, info = simplify_surface(source, 6000)
    assert info['export_allowed'], info
    good = [a for a in info['attempts'] if a['accepted']]
    assert good
    assert int(__import__('math').ceil(info['output_triangles']/500)) == min(
        int(__import__('math').ceil(a['triangles']/500)) for a in good)
    assert info['triangle_quality']['mean_ratio_p01'] >= min(.05, triangle_statistics(source)['mean_ratio_p01'])
    assert info['relative_area_change'] <= .05
    assert not info['globally_optimal']
    assert 'centroids' in info['distance_measure']
    assert 4000 <= reduced.GetNumberOfPolys() <= 8000

def test_handle_is_not_closed():
    torus = vtk.vtkParametricTorus()
    torus.SetRingRadius(70)
    torus.SetCrossSectionRadius(20)
    sampler = vtk.vtkParametricFunctionSource()
    sampler.SetParametricFunction(torus)
    sampler.SetUResolution(150)
    sampler.SetVResolution(100)
    sampler.Update()
    # Parametric source duplicates its periodic seam; clean exactly matching points.
    # Use an analytic torus sampled on a periodic grid instead.
    import numpy as np
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    u,v=np.meshgrid(np.arange(150)*2*np.pi/150,np.arange(100)*2*np.pi/100,indexing='ij')
    pts=np.stack(((70+20*np.cos(v))*np.cos(u),(70+20*np.cos(v))*np.sin(u),20*np.sin(v)),axis=-1).reshape(-1,3)
    ids=np.arange(15000).reshape(150,100)
    a=ids.ravel(); b=np.roll(ids,-1,axis=0).ravel(); c=np.roll(np.roll(ids,-1,axis=0),-1,axis=1).ravel(); d=np.roll(ids,-1,axis=1).ravel()
    faces=np.concatenate((np.stack((a,b,c),axis=1),np.stack((a,c,d),axis=1)))
    poly=vtk.vtkPolyData(); points=vtk.vtkPoints();points.SetData(numpy_to_vtk(pts,deep=True));poly.SetPoints(points)
    cells=vtk.vtkCellArray();cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0,3*len(faces)+1,3,dtype=np.int64),deep=True),numpy_to_vtkIdTypeArray(faces.ravel().astype(np.int64),deep=True));poly.SetPolys(cells)
    result,info=simplify_surface(poly,6000)
    assert info['export_allowed'], info
    assert result.GetNumberOfPoints() - result.GetNumberOfPolys()//2 == 0

def test_invalid_optimization_limits():
    with pytest.raises(ValueError):
        simplify_surface(sphere(20,20),6000,minimum_triangles=7000)
    with pytest.raises(ValueError):
        simplify_surface(sphere(20,20),max_volume_change=float('nan'))

def test_review_worker_preserves_source(tmp_path):
    import hashlib, json, subprocess
    source=sphere(120,100)
    original=tmp_path/'original_surface.vtp'
    writer=vtk.vtkXMLPolyDataWriter();writer.SetFileName(str(original));writer.SetInputData(source);assert writer.Write()
    before=hashlib.sha256(original.read_bytes()).hexdigest()
    (tmp_path/'job.json').write_text(json.dumps(dict(mode='review',target_triangles=6000,max_deviation_mm=.5)))
    worker=Path(__file__).resolve().parents[1]/'DermaDose_Surface/DermaDose_Surface/slicer_extension/DermaDoseSurface'/'surface_worker.py'
    completed=subprocess.run([sys.executable,str(worker),str(tmp_path)],capture_output=True,text=True,timeout=60)
    assert completed.returncode == 0, completed.stdout+completed.stderr
    info=json.loads((tmp_path/'result.json').read_text())
    assert info['export_allowed'],info
    assert info['max_deviation_tolerance_mm']==.5
    assert hashlib.sha256(original.read_bytes()).hexdigest()==before
    assert (tmp_path/'surface.vtp').exists()

def test_face_interior_distance_is_checked():
    import numpy as np
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    from surface_processing import shape_difference
    points=vtk.vtkPoints();points.SetData(numpy_to_vtk(np.array([[0.,0,0],[1.,0,0],[1.,1,1],[0.,1,0]]),deep=True))
    def mesh(faces):
        poly=vtk.vtkPolyData();poly.SetPoints(points)
        cells=vtk.vtkCellArray();cells.SetData(numpy_to_vtkIdTypeArray(np.array([0,3,6],dtype=np.int64),deep=True),numpy_to_vtkIdTypeArray(np.array(faces,dtype=np.int64).ravel(),deep=True));poly.SetPolys(cells)
        return poly
    # Same vertices, different non-planar diagonals: vertex-only distance is zero.
    assert shape_difference(mesh([[0,1,2],[0,2,3]]),mesh([[0,1,3],[1,2,3]])) > .1
