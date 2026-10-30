import json
from pathlib import Path
from types import SimpleNamespace
import pytest
from skindose_verifier.mrc200 import validate_configuration
from skindose_verifier.domain import InputError
from skindose_verifier.azurion import SPECIFICATION

def equipment():
    d=json.loads((Path(__file__).parents[1]/'config/equipo_mrc200_conjunto_investigacion.json').read_text())
    return SimpleNamespace(raw=d)

def test_filtration_components_not_double_counted():
    assert SPECIFICATION['tube_inherent_al_equivalent_mm']+SPECIFICATION['included_filter_plate_al_mm']==SPECIFICATION['minimum_tube_assembly_total_al_equivalent_mm']

@pytest.mark.parametrize('al',[.8,-1,float('nan')])
def test_invalid_fixed_filtration(al):
    e=equipment();e.raw['tube_characterization']['fixed_effective_al_mm']=al
    with pytest.raises(InputError):validate_configuration([],e)

def test_reject_inconsistent_table_filtration():
    e=equipment();e.raw['surface_factors']['inherent_al_mm']=3.1
    with pytest.raises(InputError):validate_configuration([],e)

def test_assembly_limit_not_insert_rating():
    with pytest.raises(InputError):validate_configuration([SimpleNamespace(kvp=150,uid='test')],equipment())
    validate_configuration([SimpleNamespace(kvp=125,uid='test')],equipment())

def test_new_dataset_provenance_and_profile_load():
    from skindose_verifier.config import Equipment
    root=Path(__file__).parents[1]
    for name in ('conjunto','conjunto_colimador'):
        Equipment.load(root/f'config/equipo_mrc200_{name}_investigacion.json')
    d=json.loads((root/'skindose_verifier/data/derived_mrc200_nominal/metadata.json').read_text())
    assert d['tables']['clinical']['rows']==1400
    assert d['clinical_validation'] is False
    assert {a['anode_angle_deg'] for a in d['spectra']}=={11}
    assert max(a['kvp'] for a in d['spectra'])==125

def test_table_angle_mismatch_is_rejected():
    from skindose_verifier.surface_factors import EventSurfaceFactors,WaterSurfaceDataset
    root=Path(__file__).parents[1]
    spec=equipment().raw['surface_factors'];spec['expected_anode_angle_deg']=12
    dataset=WaterSurfaceDataset(root/'skindose_verifier/data/derived_mrc200_nominal')
    with pytest.raises(InputError):EventSurfaceFactors(spec,dataset)
