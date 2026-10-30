from copy import deepcopy
from types import SimpleNamespace
from pathlib import Path
import json
import numpy as np
import pytest
from scipy.sparse import csr_matrix
from skindose_verifier.config import Equipment
from skindose_verifier.domain import InputError
from skindose_verifier.cli import run
from skindose_verifier.rdsr import read_rdsr
from skindose_verifier.synthetic import make_synthetic_rdsr, item
from skindose_verifier.metrics import _mesh_neighbourhood_psd
import pydicom
ROOT=Path(__file__).resolve().parents[1]

def config():
    return json.loads((ROOT/'skindose_verifier/data/equipo_sintetico.json').read_text())

@pytest.mark.parametrize('value',[-.1,float('nan'),float('inf')])
def test_rdsr_uncertainty_must_be_finite_nonnegative(value):
    d=config();d['rdsr_relative_standard_uncertainty']=value
    with pytest.raises(InputError,match='Incertidumbre'):Equipment.from_dict(d)

def test_mismatched_device_cannot_silently_replace_geometry(tmp_path):
    r=make_synthetic_rdsr(tmp_path/'r.dcm');d=config();d['model']='wrong device'
    p=tmp_path/'c.json';p.write_text(json.dumps(d))
    with pytest.raises(InputError,match='no coincide'):run(r,p,'free_air',tmp_path/'out')

def test_outside_table_cannot_silently_use_synthetic_factors(tmp_path):
    r=make_synthetic_rdsr(tmp_path/'r.dcm');ds=pydicom.dcmread(r)
    for e in ds.ContentSequence[0].ContentSequence[:2]:
        for n in e.ContentSequence:
            if n.ConceptNameCodeSequence[0].CodeValue in {'113789','113788'}:
                n.MeasuredValueSequence[0].NumericValue='20'
    ds.save_as(r,enforce_file_format=True)
    p=tmp_path/'c.json';p.write_text((ROOT/'examples/equipo_benmakhlouf_investigacion.json').read_text())
    with pytest.raises(InputError,match='fuera del dominio'):run(r,p,'skin_dose',tmp_path/'out')

def test_nonuniform_filter_not_replaced_by_arithmetic_mean(tmp_path):
    r=make_synthetic_rdsr(tmp_path/'r.dcm');ds=pydicom.dcmread(r)
    e=ds.ContentSequence[0].ContentSequence[0]
    f=next(n for n in e.ContentSequence if n.ConceptNameCodeSequence[0].CodeValue=='113771')
    f.ContentSequence.append(item('NUM','113773','maximum',value=.4,unit='mm'))
    ds.save_as(r,enforce_file_format=True)
    with pytest.raises(InputError,match='no se sustituye'):read_rdsr(r)

def test_maximum_area_average_can_have_zero_dose_center():
    # A low-dose gap may be the centre of the maximal averaging neighbourhood.
    graph=csr_matrix([[0,4,0],[4,0,4],[0,4,0]])
    s=SimpleNamespace(mesh=object(),face_graph=graph,area_mm2=np.array([33.,34.,33.]),
        points_mm=np.array([[-4,0,0],[0,0,0],[4,0,0]]))
    r=_mesh_neighbourhood_psd(s,np.array([10.,0.,10.]),100.)
    assert r['dose_gy']==pytest.approx(6.6)
    assert r['patient_point_mm']==[0,0,0]

def test_generated_dataset_declares_combined_policy(tmp_path):
    from skindose_verifier.spectral import build_dataset
    m=build_dataset(tmp_path/'table',kvps=[80],al_families=[3.1],copper=[.1],fields=[10])
    assert m['interpolation_quantity']=='combined_C'
    assert m['version']=='1.1'

def test_python_api_cannot_bypass_azurion_characterization(tmp_path):
    from skindose_verifier.engine import calculate
    d=config();d['azurion_profile']={'name':'philips_azurion'}
    e=read_rdsr(make_synthetic_rdsr(tmp_path/'r.dcm')).events
    with pytest.raises(InputError,match='caracterización'):
        calculate(e,Equipment.from_dict(d),'free_air')

def test_research_registration_does_not_mutate_caller(tmp_path):
    from skindose_verifier.engine import calculate
    d=config();d['azurion_profile']={'name':'philips_azurion',
        'allow_research_assumptions':True,'reference_table_from_first_event':True}
    eq=Equipment.from_dict(d);before=deepcopy(eq.raw)
    e=read_rdsr(make_synthetic_rdsr(tmp_path/'r.dcm')).events
    result=calculate(e,eq,'free_air')
    assert result.summary['psd_point_gy']>0
    assert eq.raw==before

def test_synthetic_quality_flag_cannot_discard_declared_aluminium():
    from skindose_verifier.surface_factors import EventSurfaceFactors,benmakhlouf_2011
    d=json.loads((ROOT/'examples/equipo_benmakhlouf_investigacion.json').read_text())
    spec=d['surface_factors'];spec['aluminium_policy']='reject_additional'
    event=SimpleNamespace(uid='test',quality={'synthetic_profile_spectrum':'additional_al_not_modeled'},
        original={'filter_stack':[{'material_scheme':'SRT','material_code':'C-120F9','selected_mm':1.}]})
    with pytest.raises(InputError,match='adicional'):
        EventSurfaceFactors(spec,benmakhlouf_2011()).evaluate_event(event,np.array([100.]))

@pytest.mark.parametrize('reference',[float('nan'),float('inf')])
def test_nonfinite_commercial_reference_is_rejected(reference):
    from skindose_verifier.validation import compare_with_radimetrics
    ours={'stage':'skin_dose','quantity':'entrance_absorbed_dose_in_water_Gy',
          'total_rdsr_ka_r_gy':1.,'psd_point_gy':1.}
    with pytest.raises(InputError):compare_with_radimetrics(ours,reference,1.,'water')
