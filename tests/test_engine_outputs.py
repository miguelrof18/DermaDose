"""Engine output contract, sparse event maps, bounded cache and configuration warnings."""
import json
from pathlib import Path

import numpy as np
import pytest

from skindose_verifier import engine
from skindose_verifier.config import Equipment
from skindose_verifier.engine import EventMaps, calculate
from skindose_verifier.rdsr import read_rdsr
from skindose_verifier.synthetic import make_synthetic_rdsr

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def setup(tmp_path):
    cfg = json.loads((ROOT/"examples/equipo_benmakhlouf_investigacion.json").read_text(encoding="utf-8"))
    events = read_rdsr(make_synthetic_rdsr(tmp_path/"synthetic.dcm", event_count=5)).events
    return cfg, events


def test_event_maps_behave_like_a_list_of_dense_maps():
    maps = EventMaps(6)
    dense = [np.array([0., 1., 0., 2., 0., 0.]), np.zeros(6), np.array([3., 0., 0., 0., 0., 4.])]
    for item in dense:
        maps.append(item)
    assert len(maps) == 3
    for stored, original in zip(maps, dense):
        np.testing.assert_array_equal(stored, original)
    np.testing.assert_array_equal(np.asarray(maps), np.stack(dense))
    np.testing.assert_array_equal(np.sum(maps, axis=0), np.sum(dense, axis=0))
    assert maps.value_at(0, 3) == 2. and maps.value_at(1, 3) == 0. and maps.value_at(2, 5) == 4.
    assert [len(m) for m in maps[1:]] == [6, 6]
    with pytest.raises(ValueError):
        maps.append(np.zeros(5))


def test_explicit_maximum_names_and_phase_timings(setup):
    cfg, events = setup
    equipment = Equipment.from_dict(cfg, base_directory=ROOT/"examples")
    for stage, name in (("free_air", "max_air_kerma_gy"), ("calibrated_air", "max_air_kerma_gy"),
                        ("skin_dose", "max_modeled_surface_dose_gy")):
        summary = calculate(events, equipment, stage, detailed=False).summary
        assert summary[name] == summary["psd_point_gy"]
        other = {"max_air_kerma_gy", "max_modeled_surface_dose_gy"} - {name}
        assert not other & summary.keys()
        assert summary["psd_point_definition"]
        for key in ("projection_seconds", "dosimetry_seconds", "metrics_seconds"):
            assert summary["performance"][key] >= 0


def test_memory_bounded_cache_does_not_change_the_map(setup, monkeypatch):
    cfg, events = setup
    equipment = Equipment.from_dict(cfg, base_directory=ROOT/"examples")
    reference = calculate(events, equipment, "skin_dose", detailed=False, cache_geometry=False)
    monkeypatch.setattr(engine, "GEOMETRY_CACHE_BYTES", 1)
    bounded = calculate(events, equipment, "skin_dose", detailed=False)
    np.testing.assert_array_equal(reference.total_gy, bounded.total_gy)
    assert bounded.summary["performance"]["geometry_cache_entries"] == 1
    np.testing.assert_array_equal(reference.total_gy, np.sum(bounded.event_maps_gy, axis=0))


def test_unknown_configuration_keys_are_reported(setup):
    cfg, events = setup
    cfg["rdsr_relative_standard_uncertainity"] = .05
    cfg["_comentario"] = "permitido"
    equipment = Equipment.from_dict(cfg, base_directory=ROOT/"examples")
    assert len(equipment.warnings) == 1 and "uncertainity" in equipment.warnings[0]
    summary = calculate(events, equipment, "free_air", detailed=False).summary
    assert summary["configuration_warnings"] == list(equipment.warnings)
