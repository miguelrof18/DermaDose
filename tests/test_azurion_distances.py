from dataclasses import replace
from types import SimpleNamespace
import pytest

from skindose_verifier.azurion import prepare_events
from skindose_verifier.domain import InputError, IrradiationEvent


def equipment(policy='reject'):
    geometry = dict(source_to_detector_mm=1044.3, source_to_isocenter_mm=810., source_to_irp_mm=660.)
    return SimpleNamespace(raw={'azurion_profile': dict(name='philips_azurion',
        allow_research_assumptions=True, invalid_distance_policy=policy)},
        geometry_for=lambda _: geometry)


def event(**kwargs):
    return replace(IrradiationEvent(uid='test', plane='single', sid_mm=996.7,
        sod_mm=810., srp_mm=660., ka_r_gy=.00480598), **kwargs)


def test_strict_error_identifies_event_and_keeps_nonzero_dose():
    with pytest.raises(InputError, match=r'evento 2/2 .*sid_mm=0.*sod_mm=0.*srp_mm=0.*0.00480598'):
        prepare_events([event(), event(sid_mm=0, sod_mm=0, srp_mm=0)], equipment())


def test_research_changes_only_invalid_distances_and_audits_originals():
    original = event(sid_mm=0, sod_mm=None, srp_mm=0)
    valid = event()
    results = prepare_events([valid, original], equipment('configured_research'))
    assert results[0].sid_mm == 996.7
    assert results[1].sid_mm == 1044.3
    assert results[1].sod_mm == 810.
    assert results[1].srp_mm == 660.
    assert results[1].ka_r_gy == original.ka_r_gy
    assert results[1].original['research_distance_replacements']['sid_mm']['original_mm'] == 0
    assert results[1].original['research_distance_replacements']['sod_mm']['original_mm'] is None
    assert original.sid_mm == 0
    assert not original.quality


def test_research_rejects_invalid_configured_distance():
    eq = equipment('configured_research')
    eq.geometry_for = lambda _: {'source_to_detector_mm': 0}
    with pytest.raises(InputError, match='JSON'):
        prepare_events([event(sid_mm=0)], eq)
