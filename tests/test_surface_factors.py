"""Published nodes, bounded interpolation, audit and full dose-chain integration."""
from copy import deepcopy
import csv
from dataclasses import replace
from importlib.resources import files
import json
from pathlib import Path

import numpy as np
import pydicom
import pytest

from skindose_verifier.cli import run
from skindose_verifier.config import Equipment
from skindose_verifier.domain import InputError
from skindose_verifier.engine import calculate
from skindose_verifier.factor_tables import RectilinearFactorTable
from skindose_verifier.geometry import project
from skindose_verifier.patient import sample_patient
from skindose_verifier.rdsr import read_rdsr
from skindose_verifier.surface_factors import benmakhlouf_2011, WaterSurfaceDataset
from skindose_verifier.synthetic import make_synthetic_rdsr

ROOT = Path(__file__).resolve().parents[1]
DATA = files('skindose_verifier').joinpath('data/benmakhlouf_2011')


@pytest.fixture
def setup(tmp_path):
    cfg=json.loads((ROOT/'examples/equipo_benmakhlouf_investigacion.json').read_text())
    events=read_rdsr(make_synthetic_rdsr(tmp_path/'synthetic.dcm')).events
    return cfg,events


@pytest.mark.parametrize('kind,file', [('monoenergetic','table1_monoenergetic.csv'),('clinical','table3_clinical.csv')])
def test_every_published_pair_is_recovered_exactly(kind,file):
    dataset=benmakhlouf_2011()
    rows=list(csv.DictReader(DATA.joinpath(file).read_text().splitlines()))
    for row in rows:
        args=({k:float(row[k]) for k in ('kvp','added_cu_mm','inherent_al_mm','field_side_cm')} if kind=='clinical'
              else {k:float(row[k]) for k in ('energy_kev','field_side_cm')})
        result=dataset.evaluate(kind,**args)
        assert float(result.k_bs)==float(row['k_bs'])
        assert float(result.k_f)==float(row['k_f'])
        assert result.audit['status']=='tabulated'
        assert result.audit['source_nodes']==[row]
        assert result.audit['extrapolated'] is False


def test_bilinear_monoenergetic_values_by_hand():
    r=benmakhlouf_2011().evaluate('monoenergetic',energy_kev=42.5,field_side_cm=15)
    assert float(r.k_bs)==pytest.approx((1.403+1.469+1.447+1.543)/4,abs=1e-15)
    assert float(r.k_f)==pytest.approx((1.016+1.016+1.020+1.020)/4,abs=1e-15)
    assert len(r.audit['source_nodes'])==4


def test_trilinear_clinical_values_by_hand():
    r=benmakhlouf_2011().evaluate('clinical',inherent_al_mm=3.1,kvp=85,added_cu_mm=.15,field_side_cm=15)
    assert float(r.k_bs)==pytest.approx(sum([1.392,1.476,1.421,1.522,1.401,1.495,1.426,1.536])/8,abs=1e-15)
    assert float(r.k_f)==pytest.approx(sum([1.025,1.024,1.028,1.027,1.029,1.028,1.033,1.032])/8,abs=1e-15)
    assert len(r.audit['source_nodes'])==8
    assert r.audit['status']=='interpolated'


@pytest.mark.parametrize('side,expected',[(10,1.392),(12.5,1.413),(20,1.476)])
def test_field_interpolation_is_in_side_not_area(side,expected):
    r=benmakhlouf_2011().evaluate('clinical',inherent_al_mm=3.1,kvp=80,added_cu_mm=.1,field_side_cm=side)
    assert float(r.k_bs)==pytest.approx(expected)


@pytest.mark.parametrize('axis,value',[
 ('kvp',49.999),('kvp',150.001),('added_cu_mm',-.001),('added_cu_mm',.901),
 ('field_side_cm',4.999),('field_side_cm',35.001),('kvp',float('nan')),('field_side_cm',float('inf'))])
def test_outside_domain_never_extrapolates(axis,value):
    args=dict(inherent_al_mm=3.1,kvp=80,added_cu_mm=.1,field_side_cm=20)
    args[axis]=value
    with pytest.raises(InputError,match='no se extrapola'):
        benmakhlouf_2011().evaluate('clinical',**args)


def test_family_specific_domain_and_no_invented_al_or_low_energy():
    d=benmakhlouf_2011()
    for al,cu in [(4.1,.3),(5.1,.001),(3.0,.1)]:
        with pytest.raises(InputError):
            d.evaluate('clinical',inherent_al_mm=al,kvp=80,added_cu_mm=cu,field_side_cm=20)
    with pytest.raises(InputError):
        d.evaluate('monoenergetic',energy_kev=4,field_side_cm=10)
    assert d.evaluate('clinical',inherent_al_mm=5.1,kvp=80,added_cu_mm=0,field_side_cm=20).k_bs==1.446


def test_bad_h_preserved_and_flagged_but_not_used_as_axis():
    r=benmakhlouf_2011().evaluate('clinical',inherent_al_mm=3.1,kvp=120,added_cu_mm=0,field_side_cm=5)
    assert r.k_bs==1.255 and r.k_f==1.034
    assert r.audit['source_nodes'][0]['h_original']=='3.29'
    assert r.audit['source_quality_flags']==['h_outside_0_1_in_source']


@pytest.mark.parametrize('key,value',[
 ('target','Mo'),('inherent_al_mm',None),('out_of_domain','clip'),('incident_kerma_basis','surface_air_kerma'),
 ('allow_water_phantom_approximation',False),('quality_basis','hvl_only'),('field_shape_policy','guess')])
def test_invalid_provider_configuration_fails(setup,key,value):
    cfg,_=setup;cfg['surface_factors'][key]=value
    with pytest.raises(InputError): Equipment.from_dict(cfg)


def test_water_medium_and_no_double_counting_including_modes(setup):
    cfg,_=setup
    wrong=deepcopy(cfg);wrong['dose_medium']='skin'
    with pytest.raises(InputError,match='water'): Equipment.from_dict(wrong)
    for name in ('bsf','f_skin'):
        for mode in (False,True):
            wrong=deepcopy(cfg)
            block=wrong.setdefault('factors_by_mode',{}).setdefault('fluoro',{}) if mode else wrong['factors']
            block[name]={'value':1.1,'identifier':'legacy','provenance':'test'}
            with pytest.raises(InputError,match='doble conteo'): Equipment.from_dict(wrong)


def test_added_aluminium_from_rdsr_is_not_silently_folded_into_inherent_al(setup):
    cfg, events = setup
    event = replace(events[0], original={**events[0].original, 'filter_stack': [
        {'material_scheme': 'SRT', 'material_code': 'C-120F9', 'selected_mm': 1.0}
    ]})
    provider = Equipment.from_dict(cfg).surface_factor_provider
    with pytest.raises(InputError, match='filtración adicional de Al'):
        provider.evaluate_event(event, np.asarray([100.0]))


def test_out_of_domain_uses_whole_run_historical_synthetic_factors(tmp_path):
    config = tmp_path / 'benmakhlouf.json'
    data=json.loads((ROOT/'examples/equipo_benmakhlouf_investigacion.json').read_text(encoding='utf-8'))
    data['allow_historical_factor_fallback']=True
    config.write_text(json.dumps(data),encoding='utf-8')
    rdsr = make_synthetic_rdsr(tmp_path / 'small_field.dcm')
    ds = pydicom.dcmread(rdsr)
    for event in ds.ContentSequence[0].ContentSequence[:2]:
        for node in event.ContentSequence:
            if node.ConceptNameCodeSequence[0].CodeValue in {'113789', '113788'}:
                node.MeasuredValueSequence[0].NumericValue = '20'
    ds.save_as(rdsr, enforce_file_format=True)
    result = run(rdsr, config, 'skin_dose', tmp_path / 'result')
    summary = json.loads(Path(result['files']['summary']).read_text(encoding='utf-8'))
    assert summary['status'] == 'research_synthetic_factors_out_of_domain_no_clinical_use'
    assert summary['audit']['effective_factor_source'] == 'historical_synthetic_whole_run'
    assert summary['audit']['factor_fallback_reason'] is not None
    assert summary['factor_metadata']['default']['bsf']['value'] == 1.25
    assert summary['factor_metadata']['default']['f_skin']['value'] == 1.05


@pytest.mark.parametrize('kvp,cu,side,bs,f',[(80,.1,20,1.476,1.024),(50,0,5,1.204,1.016),
 (150,.9,35,1.592,1.067),(85,.15,15,1.458625,1.02825)])
def test_final_dose_contains_exactly_one_pair_with_geometry_and_table(setup,kvp,cu,side,bs,f):
    cfg,events=setup;eq=Equipment.from_dict(cfg)
    # Flat surface is 600 mm from focus, detector at 1100 mm.
    width=side*10*1100/600
    ev=replace(events[0],kvp=kvp,filtration_mm_cu=cu,field_width_mm=width,field_height_mm=width)
    air=calculate([ev],eq,'calibrated_air',detailed=False)
    dose=calculate([ev],eq,'skin_dose',detailed=False)
    np.testing.assert_allclose(dose.total_gy,air.total_gy*bs*f,rtol=2e-15)
    mask=dose.total_gy>0
    distance2=np.sum((dose.surface.points_mm[mask]-[0,0,-700])**2,axis=1)
    expected=ev.ka_r_gy*550**2/distance2*.8*bs*f
    np.testing.assert_allclose(dose.total_gy[mask],expected,rtol=2e-15)
    audit=dose.summary['events'][0]['surface_factor_audit_per_sample'][0]
    assert audit['table']==3 and audit['clinical_validation'] is False
    assert dose.summary['quantity']=='entrance_absorbed_dose_in_water_Gy'
    assert 'incertidumbre de bsf' in dose.summary['uncertainty']['unquantified_components']
    json.dumps(dose.summary,allow_nan=False)


def test_per_event_factors_and_accumulation(setup):
    cfg,events=setup;eq=Equipment.from_dict(cfg)
    a=replace(events[0],kvp=80,filtration_mm_cu=.1)
    b=replace(a,uid='second',kvp=100,filtration_mm_cu=.3)
    both=calculate([a,b],eq,'skin_dose',detailed=False)
    one=calculate([a],eq,'skin_dose',detailed=False)
    two=calculate([b],eq,'skin_dose',detailed=False)
    np.testing.assert_array_equal(both.total_gy,one.total_gy+two.total_gy)
    assert both.summary['events'][0]['factor_ranges_per_sample']!=both.summary['events'][1]['factor_ranges_per_sample']


def test_rectangular_requires_explicit_mapping_and_errors_name_event(setup):
    cfg,events=setup;ev=replace(events[0],field_height_mm=300)
    with pytest.raises(InputError,match='campo no cuadrado'):
        calculate([ev],Equipment.from_dict(cfg),'skin_dose')
    cfg['surface_factors']['field_shape_policy']='equal_area'
    r=calculate([ev],Equipment.from_dict(cfg),'skin_dose',detailed=False)
    assert r.summary['events'][0]['surface_factor_audit_per_sample'][0]['field_mapping'].startswith('equal_area')
    bad=replace(ev,uid='OUTSIDE-123',kvp=151)
    with pytest.raises(InputError,match='OUTSIDE-123'):
        calculate([bad],Equipment.from_dict(cfg),'skin_dose')


def test_custom_dataset_separate_and_checksum_enforced(setup,tmp_path):
    import shutil
    target=tmp_path/'own'
    shutil.copytree(str(DATA),target)
    m=json.loads((target/'metadata.json').read_text());m['dataset_id']='MY-INDEPENDENT-DATA';m['version']='test'
    (target/'metadata.json').write_text(json.dumps(m))
    cfg,events=setup
    cfg['surface_factors'].update(provider='tabulated_water',dataset_directory='own')
    eq=Equipment.from_dict(cfg,base_directory=tmp_path)
    assert eq.surface_factor_provider.dataset.metadata['dataset_id']=='MY-INDEPENDENT-DATA'
    assert benmakhlouf_2011().metadata['dataset_id']=='KS-ASF-201101-IR'
    calculate([events[0]],eq,'skin_dose',detailed=False)
    file=target/'table3_clinical.csv';file.write_text(file.read_text().replace('1.192','1.193',1))
    with pytest.raises(InputError,match='SHA-256'): WaterSurfaceDataset(target)


def test_generic_interpolator_rejects_missing_duplicate_and_invalid_nodes():
    rows=[dict(row_id=str(i),x=x,y=y,k_bs=1+x+y,k_f=1) for i,(x,y) in enumerate([(0,0),(1,0),(0,1),(1,1)])]
    for invalid in (rows[:-1],rows+[rows[0]],[{**r,'k_bs':float('nan')} for r in rows]):
        with pytest.raises(InputError): RectilinearFactorTable(invalid,['x','y'])
    grid=RectilinearFactorTable(rows,['x','y'])
    r=grid.evaluate(x=np.array([0,.5,1]),y=.5)
    np.testing.assert_array_equal(r.k_bs,[1.5,2,2.5])
    with pytest.raises(InputError): grid.evaluate(x=.5)
