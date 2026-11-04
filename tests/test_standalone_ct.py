import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
import vtk
ROOT=Path(__file__).resolve().parents[1]/'DermaDose_Surface'/'standalone_ct'
sys.path.insert(0,str(ROOT))
from dermadose_ct import segment_body, make_surface
from surface_processing import surface_quality


def phantom():
    z,y,x=np.indices((60,70,80))
    body=((x-40)/22)**2+((y-34)/23)**2+((z-30)/25)**2<1
    data=np.full(body.shape,-1000,dtype=np.int16);data[body]=40
    data[20:39,28:35,33:48]=-1000
    data[3:6,3:6,3:6]=80
    data[:,64:67,:]=0  # separate table
    return data


def test_body_holes_table_and_no_input_mutation():
    data=phantom();original=data.copy()
    mask,audit=segment_body(data)
    assert np.array_equal(data,original)
    assert mask[30,30,40]==1
    assert not mask[:,64:67].any()
    assert not mask[4,4,4]
    assert audit['discarded_voxels']>0
    assert audit['filled_voxels']>0


def test_connected_table_not_claimed_removed():
    data=phantom();data[28:33,55:65,38:43]=40
    mask,audit=segment_body(data)
    assert mask[:,64:67].any()
    assert any('mesa' in w for w in audit['warnings'])


def test_hollow_airway_and_slice_option():
    data=np.full((24,24,24),-1000,dtype=np.int16)
    data[3:21,3:21,3:21]=40;data[:,10:14,10:14]=-1000
    without,_=segment_body(data,fill_slices=False)
    with_fill,_=segment_body(data,fill_slices=True)
    assert not without[12,12,12] and with_fill[12,12,12]


def test_oblique_affine_crop_offset_exact():
    mask=np.zeros((30,40,50),dtype=np.uint8);mask[10:20,15:30,20:40]=1
    matrix=np.array([[0,-2,0,100],[1.5,0,0,-40],[0,0,3,200],[0,0,0,1.]])
    poly=make_surface(mask,matrix)
    expected=[100-2*29.5,100-2*14.5,-40+1.5*19.5,-40+1.5*39.5,200+3*9.5,200+3*19.5]
    assert np.allclose(poly.GetBounds(),expected)
    quality=surface_quality(poly)
    assert quality['boundary_edges']==quality['nonmanifold_edges']==0


def test_image_edge_cap_reported():
    data=np.full((20,20,20),-1000,dtype=np.int16);data[:15,4:16,4:16]=50
    mask,audit=segment_body(data)
    assert 'K−' in audit['touches_image_faces']
    assert surface_quality(make_surface(mask,np.eye(4)))['boundary_edges']==0


def test_bad_input_rejected():
    with pytest.raises(ValueError):segment_body(np.zeros((5,5)))
    with pytest.raises(ValueError):segment_body(np.full((5,5,5),-1000))
    with pytest.raises(ValueError):make_surface(np.ones((5,5,5)),np.zeros((4,4)))


def test_standalone_worker_segment_and_mesh(tmp_path):
    np.save(tmp_path/'input.npy',phantom())
    job=dict(action='segment',threshold_hu=-350,fill_slices=True,ijk_to_ras_mm=np.eye(4).tolist())
    def run():
        (tmp_path/'job.json').write_text(json.dumps(job))
        p=subprocess.run([sys.executable,str(ROOT/'dermadose_ct.py'),'--worker',str(tmp_path)],capture_output=True,text=True,timeout=90)
        assert p.returncode==0,p.stdout+p.stderr
    run()
    np.save(tmp_path/'input.npy',np.load(tmp_path/'mask.npy'))
    job['action']='mesh';run()
    result=json.loads((tmp_path/'result.json').read_text())
    assert result['export_allowed'],result
    assert result['output_triangles']<=200000
    assert (tmp_path/'surface.stl').exists()
    assert (tmp_path/'preview.vtp').exists()
    assert not result['registration_established']
