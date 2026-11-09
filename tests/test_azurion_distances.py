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


def test_old_research_policy_no_longer_replaces_distances():
    original = event(sid_mm=0, sod_mm=None, srp_mm=0)
    with pytest.raises(InputError, match='No hay una recuperación admisible'):
        prepare_events([original], equipment('configured_research'))
    assert original.sid_mm == 0
    assert not original.quality


def test_valid_distances_are_preserved():
    original = event()
    result, = prepare_events([original], equipment())
    assert result.sid_mm == original.sid_mm
    assert result.sod_mm == original.sod_mm
    assert result.srp_mm == original.srp_mm


def test_removed_research_option_is_rejected():
    from skindose_verifier.example_runner import main
    with pytest.raises(SystemExit) as exc:
        main(['--research-distances'])
    assert exc.value.code == 2
