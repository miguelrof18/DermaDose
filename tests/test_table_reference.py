from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import json
import numpy as np

from skindose_verifier.config import Equipment
from skindose_verifier.engine import normalize, calculate
from skindose_verifier.geometry import table_shift
from skindose_verifier.rdsr import read_rdsr
from skindose_verifier.synthetic import make_synthetic_rdsr


def test_isocenter_to_legacy_preserves_selected_axes_and_physical_shift(tmp_path):
    base = Equipment.load(Path(__file__).parents[1]/'skindose_verifier/data/equipo_sintetico.json')
    raw = deepcopy(base.raw)
    raw['geometry'].update(table_coordinates='legacy', table_origin_mm=[300,400,500],
                           table_matrix=[[0,1,0],[-1,0,0],[0,0,1]], isocenter_mm=[10,20,30])
    equipment = Equipment.from_dict(raw)
    event = read_rdsr(make_synthetic_rdsr(tmp_path/'source.dcm')).events[0]
    event = replace(event, table_mm=(5,6,7))
    converted = normalize(event, equipment)
    assert converted.table_system == 'legacy'
    np.testing.assert_allclose(table_shift(converted, equipment), [16,15,37])
    assert converted.original['table_coordinate_conversion']['source_mm'] == [5,6,7]
    assert event.table_system == 'isocenter'
    assert event.table_mm == (5,6,7)
    # A second normalization must not repeat the translation.
    assert normalize(converted, equipment).table_mm == converted.table_mm


def test_equivalent_references_give_identical_dose(tmp_path):
    base = Equipment.load(Path(__file__).parents[1]/'skindose_verifier/data/equipo_sintetico.json')
    events = read_rdsr(make_synthetic_rdsr(tmp_path/'source.dcm')).events
    raw = deepcopy(base.raw)
    raw['geometry'].update(table_coordinates='legacy', table_origin_mm=[300,400,500])
    adapted = Equipment.from_dict(raw)
    expected = calculate(events, base)
    actual = calculate(events, adapted)
    np.testing.assert_allclose(actual.total_gy, expected.total_gy)
    assert actual.summary['events'][0]['table_coordinate_conversion']['source_system'] == 'isocenter'


def test_explicit_isocenter_axes_override_legacy_axes(tmp_path):
    base = Equipment.load(Path(__file__).parents[1]/'skindose_verifier/data/equipo_sintetico.json')
    raw = deepcopy(base.raw)
    raw['geometry'].update(table_coordinates='legacy', table_origin_mm=[10,20,30],
                           isocenter_table_matrix=[[1,0,0],[0,0,1],[0,-1,0]])
    equipment = Equipment.from_dict(raw)
    event = read_rdsr(make_synthetic_rdsr(tmp_path/'source.dcm')).events[0]
    converted = normalize(replace(event, table_mm=(10,20,30)), equipment)
    np.testing.assert_allclose(table_shift(converted, equipment), [10,30,-20])
