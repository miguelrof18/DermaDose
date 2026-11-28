"""Clear errors for stale or incoherent inputs in a case folder."""
import json

import numpy as np
import pytest

from skindose_verifier.domain import InputError
from skindose_verifier.example_runner import ExampleError, prepare_example
from skindose_verifier.placement import select_target


def _box(z_low, z_high):
    return np.array([[x, y, z] for x in (-150, 150) for y in (-100, 100)
                     for z in np.arange(z_low, z_high + 1, 5.)], float)


def test_shoulder_anchor_at_caudal_end_of_mesh_is_explained():
    vertices = _box(846., 1297.)
    registration = {"anchor_ct_lps_mm": [0., 100., 847.], "anchor_table_mm": [0., -39., -40.],
                    "rotation_table_from_ct": np.eye(3).tolist(),
                    "orientation_basis": "operator_confirmed_aligned_head_first_supine"}
    with pytest.raises(InputError, match="no es coherente con la malla"):
        select_target(vertices, {}, (), registration)
    registration["anchor_ct_lps_mm"] = [0., 100., 1290.]
    target, info = select_target(vertices, {}, (), registration)
    assert info["target_source"].startswith("cardiac_level_proxy") and 846 < target[2] < 1290


def test_pre_0_7_azurion_profile_in_case_folder_is_rejected(tmp_path):
    old = {"version": "v", "manufacturer": "Philips", "model": "Azurion", "factors": {},
           "geometry": {"table_matrix": np.eye(3).tolist(), "table_origin_mm": [0, 0, 0]},
           "azurion_profile": {"reference_table_from_first_event": True}}
    (tmp_path / "viejo.json").write_text(json.dumps(old), encoding="utf-8")
    (tmp_path / "caso.dcm").write_bytes(b"")
    with pytest.raises(ExampleError, match="anterior a DermaDose 0.7.0"):
        prepare_example(tmp_path, mesh="xcat")


def test_kerma_weighted_reference_is_the_working_position():
    from types import SimpleNamespace
    from skindose_verifier.placement import choose_reference_event

    def event(long_, lat, height, ka):
        reading = lambda v: [{"value": v}]
        return SimpleNamespace(original={"table_long": reading(long_), "table_lat": reading(lat),
                                         "table_height": reading(height)},
                               table_system="isocenter", table_mm=None, ka_r_gy=ka, uid=str(ka))
    # Low-dose scout 85 mm away, then the working position with most of the kerma.
    events = [event(-13, 859, 843, .0001), event(-58, 941, 882, .02), event(-64, 942, 882, .03),
              event(-40, 909, 862, .004)]
    assert choose_reference_event(events, "first")[0] == 1
    assert choose_reference_event(events, "kerma_weighted")[0] == 3
    assert choose_reference_event(events)[0] == 3


def test_shoulder_anchor_a_few_centimetres_from_the_head_edge_is_rejected():
    from skindose_verifier.reference_frames import check_shoulder_anchor_distance
    registration = {"anchor_table_mm": [0., -39., -40.], "anchor_definition": "posterior_shoulder_midline"}
    with pytest.raises(InputError, match="40 cm son 400 mm"):
        check_shoulder_anchor_distance(registration)
    registration["anchor_table_mm"] = [0., -39., -400.]
    check_shoulder_anchor_distance(registration)
