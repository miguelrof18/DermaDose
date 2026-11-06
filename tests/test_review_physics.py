"""Independent regression cases for the November physics review."""
from dataclasses import replace
from pathlib import Path
import numpy as np
import pytest
from skindose_verifier.config import Equipment
from skindose_verifier.domain import InputError
from skindose_verifier.engine import calculate
from skindose_verifier.geometry import beam_from_event, project
from skindose_verifier.irp import verify_irp
from skindose_verifier.patient import sample_patient
from skindose_verifier.rdsr import read_rdsr, parse_event
from skindose_verifier.synthetic import make_synthetic_rdsr, item, coded
import pydicom

ROOT = Path(__file__).resolve().parents[1]

@pytest.fixture
def case(tmp_path):
    path = make_synthetic_rdsr(tmp_path/'synthetic.dcm')
    cfg = Equipment.load(ROOT/'skindose_verifier/data/equipo_sintetico.json')
    return path, cfg, read_rdsr(path).events


def test_dicom_detector_plane_survives_isocenter_configuration(case):
    _, cfg, events = case
    ev = events[0]
    expected = project(ev, cfg, sample_patient(cfg.patient))
    expected_irp = verify_irp(ev, cfg)['kap_area']
    cfg.geometry['field_reference'] = 'isocenter'
    assert beam_from_event(ev, cfg).reference_distance_mm == ev.sid_mm
    actual = project(ev, cfg, sample_patient(cfg.patient))
    np.testing.assert_array_equal(actual.visible_mask, expected.visible_mask)
    np.testing.assert_allclose(actual.field_area_cm2, expected.field_area_cm2, rtol=0, atol=0)
    assert verify_irp(ev, cfg)['kap_area'] == expected_irp


def test_manual_event_still_supports_explicit_isocenter(case):
    _, cfg, events = case
    ev = replace(events[0], field_reference='isocenter')
    assert beam_from_event(ev, cfg).reference_distance_mm == ev.sod_mm
    cfg.geometry['field_reference'] = 'isocenter'
    assert beam_from_event(replace(ev, field_reference=None), cfg).reference_distance_mm == ev.sod_mm


def test_mixed_reference_cache_does_not_reuse_wrong_projection(case):
    _, cfg, events = case
    ev = events[0]
    pair = [ev, replace(ev, uid='another', field_reference='isocenter')]
    cached = calculate(pair, cfg, detailed=False)
    uncached = calculate(pair, cfg, detailed=False, cache_geometry=False)
    assert cached.summary['performance']['projection_cache_hits'] == 0
    np.testing.assert_array_equal(cached.total_gy, uncached.total_gy)


def test_duplicate_python_event_cannot_double_dose(case):
    _, cfg, events = case
    with pytest.raises(InputError, match='UID.*repetido'):
        calculate([events[0], events[0]], cfg, detailed=False)


@pytest.mark.parametrize('scheme,code', [('SCT','12503006'),('SRT','C-12000'),('SRT','C-120F9')])
def test_aluminium_codes_flow_through_parser_and_dosimetry(case, scheme, code):
    path, _, _ = case
    node = pydicom.dcmread(path).ContentSequence[0].ContentSequence[0]
    container = item('CONTAINER', '113771', 'X-Ray Filters')
    container.ContentSequence.extend([
        item('CODE', '113757', 'material', value_code=coded(code, 'label ignored', scheme)),
        item('NUM', '113758', 'thickness', value=1, unit='mm')])
    node.ContentSequence.append(container)
    event = parse_event(node)
    assert event.filtration_mm_cu == .2
    cfg = Equipment.load(ROOT/'examples/equipo_benmakhlouf_investigacion.json')
    provider = cfg.surface_factor_provider
    with pytest.raises(InputError, match='filtración adicional de Al'):
        provider.evaluate_event(event, np.array([100.]))
    provider.aluminium_policy = 'fixed_plus_event'
    event = replace(event, filtration_mm_cu=0.)
    result = provider.evaluate_event(event, np.array([100.]))
    assert result.audit['aluminium_mapping']['event_added_al_mm'] == 1.
    assert result.audit['aluminium_mapping']['selected_table_al_mm'] == pytest.approx(provider.al+1)


def test_negative_collimated_area_is_controlled_error(case):
    path, _, _ = case
    node = pydicom.dcmread(path).ContentSequence[0].ContentSequence[0]
    node.ContentSequence = [n for n in node.ContentSequence if n.ConceptNameCodeSequence[0].CodeValue not in ('113788','113789')]
    node.ContentSequence.append(item('NUM', '113790', 'area', value=-1, unit='m2'))
    with pytest.raises(InputError, match='área.*no positiva'):
        parse_event(node)


def test_accumulation_matches_original_order_and_retains_event_maps(case):
    _, cfg, events = case
    seq = [replace(events[i%2], uid=f'event-{i}', ka_r_gy=10.**(-i%9-6)) for i in range(31)]
    result = calculate(seq, cfg, detailed=False)
    np.testing.assert_array_equal(result.total_gy, np.sum(np.stack(result.event_maps_gy), axis=0))
    assert len(result.event_maps_gy) == 31


def test_unresolved_irradiated_neighbourhood_is_not_zero(case):
    from skindose_verifier.metrics import neighbourhood_psd
    _, cfg, _ = case
    cfg.patient['resolution_mm'] = 2.5
    surface = sample_patient(cfg.patient)
    dose = np.zeros(len(surface.points_mm))
    dose[np.linalg.norm(surface.points_mm[:, :2], axis=1) < 30] = 1.
    # Interior lattice discs cover 131.25 mm2 and fail the 125 mm2 upper limit;
    # admissible boundary discs all carry zero dose. Zero is not a peak estimate.
    assert neighbourhood_psd(surface, dose) is None
    result = neighbourhood_psd(surface, np.zeros_like(dose))
    assert result is not None and result['dose_gy'] == 0.


def test_unresolved_neighbourhood_report_status(case):
    _, cfg, events = case
    cfg.patient['resolution_mm'] = 2.5
    result = calculate(events, cfg)
    assert result.summary['psd_point_gy'] > 0
    assert result.summary['psd_neighbourhood_1cm2'] is None
    assert result.summary['psd_neighbourhood_1cm2_status'] == 'unavailable_spatial_sampling'
