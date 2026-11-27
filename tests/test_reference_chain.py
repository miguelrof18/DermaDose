"""C -> T -> I chain: the worked example of docs/sistemas_de_referencia.md,
estimated links and table calibration. NumPy only."""
from dataclasses import replace

import numpy as np
import pytest

from skindose_verifier.domain import InputError, IrradiationEvent
from skindose_verifier.placement import (choose_reference_event, estimate_ct_registration, estimate_room_reference,
                                         mesh_vertices_lps, pose_rotation, select_target)
from skindose_verifier.reference_frames import (RegisteredFrames, ct_to_table, motion_matrix, rigid,
                                                room_reference_state)
from skindose_verifier.table_calibration import calibrated_profile, parse_move

# Worked example of the reference-frame document (illustrative, not Azurion data).
A_C, A_T, P_C = [20, 110, -140], [0, -35, -280], [20, 110, -360]
S = np.array([[0, 1, 0], [0, 0, -1], [-1, 0, 0]], dtype=float)
M0, M1 = [420, -10, 155], [320, 35, 185]
Q_T, Q_I0 = [95, -5, -225], [95, 120, 275]


def event(readings, uid="e", ka=1., pose="HFS"):
    long, lat, height = readings
    original = {k: [{"value": float(v)}] for k, v in (("table_long", long), ("table_lat", lat), ("table_height", height))}
    original["intervention_patient_position"] = pose
    return IrradiationEvent(uid=uid, plane="113622", ka_r_gy=ka, table_mm=(lat, long, height),
                            table_system="legacy", original=original)


def frames(s=S, origin=None):
    transform, _ = ct_to_table({"anchor_ct_lps_mm": A_C, "anchor_table_mm": A_T,
                                "rotation_table_from_ct": np.eye(3).tolist()})
    origin = np.subtract(Q_I0, Q_T) if origin is None else np.asarray(origin, dtype=float)
    return RegisteredFrames(transform, np.asarray(M0, float), origin, np.asarray(s, float), {})


def test_document_worked_example():
    f = frames()
    np.testing.assert_allclose(f.table_from_ct[:3, 3], [-20, -145, -140])            # t_T<-C
    np.testing.assert_allclose(f.reference_origin_i, [0, 125, 500])                  # t_0
    np.testing.assert_allclose(S @ np.subtract(M1, M0), [45, -30, 100])               # S Δm_1
    for readings, origin, point in ((M0, [0, 125, 500], [0, 90, 0]), (M1, [45, 95, 600], [45, 60, 100])):
        t_k = f.table_origin_i(event(readings))
        np.testing.assert_allclose(t_k, origin)
        m_i_t = rigid(np.eye(3), t_k)
        m_i_c = m_i_t @ f.table_from_ct
        np.testing.assert_allclose(m_i_c @ [*P_C, 1], [*point, 1])
        np.testing.assert_allclose(np.linalg.inv(m_i_c) @ [*point, 1], [*P_C, 1])
    np.testing.assert_allclose((rigid(np.eye(3), [45, 95, 600]) @ f.table_from_ct)[:3, 3], [25, -50, 460])


def test_unknown_motion_column_only_blocks_when_that_reading_moves():
    partial = S.copy()
    partial[:, 1] = np.nan
    f = frames(partial)
    np.testing.assert_allclose(f.table_origin_i(event([500, -10, 155])), [0, 125, 420])
    with pytest.raises(InputError, match="lateral .* cambia \\+45.0 mm"):
        f.table_origin_i(event([420, 35, 155]))
    matrix = motion_matrix({"motion_matrix_i_per_mm": [[1, 0, 0], [0, 0, -1], [0, None, 0]]})
    assert np.isnan(matrix[2, 1]) and matrix[0, 0] == 1
    with pytest.raises(InputError, match="independientes"):
        motion_matrix({"motion_matrix_i_per_mm": [[1, 1, None], [0, 0, None], [0, 0, None]]})


def test_room_reference_states():
    assert room_reference_state({"reference": None}) == "missing"
    assert room_reference_state({"reference": {"marker_table_mm": None}}) == "missing"
    assert room_reference_state({"reference": {"readings_long_lat_height_mm": M0, "marker_table_mm": Q_T,
                                                "marker_isocenter_mm": Q_I0}}) == "measured"
    estimated = {"readings_long_lat_height_mm": M0, "origin_isocenter_mm": [0, 125, 500], "status": "estimated_x"}
    assert room_reference_state({"reference": estimated}) == "estimated"
    with pytest.raises(InputError, match="status"):
        room_reference_state({"reference": {**estimated, "status": "measured"}})
    with pytest.raises(InputError, match="no ambas"):
        room_reference_state({"reference": {**estimated, "marker_table_mm": Q_T}})


def cylinder(a=170., b=110., z0=-600., z1=0., n=64, levels=31, centre=(20., 0., 0.)):
    theta = np.linspace(-np.pi, np.pi, n, endpoint=False)
    return np.array([[a*np.cos(t), b*np.sin(t), z] for z in np.linspace(z0, z1, levels) for t in theta]) + centre


def test_estimated_registration_rests_the_back_on_the_mattress():
    vertices = cylinder()
    target = np.array([20., -10., -300.])
    rotation, pose = pose_rotation([event(M0)])
    assert pose["position"] == "HFS" and np.allclose(rotation, np.eye(3))
    registration = estimate_ct_registration(vertices, target, rotation, mattress_mm=39, target_distance_from_head_mm=580,
                                            mesh_sha256="x", target_info={}, pose_info=pose)
    np.testing.assert_allclose(registration["anchor_ct_lps_mm"], [20, 110, -300])
    transform, _ = ct_to_table(registration)
    p_t = vertices @ transform[:3, :3].T + transform[:3, 3]
    slab = np.abs(p_t[:, 2] + 580) <= 50
    assert p_t[slab, 1].max() == pytest.approx(-39)                 # back on the mattress, T y points down
    np.testing.assert_allclose(transform[:3, :3] @ target + transform[:3, 3], [0, -159, -580])


def test_prone_pose_rests_the_anterior_surface():
    vertices = cylinder()
    rotation, pose = pose_rotation([event(M0, pose="HFP")])
    np.testing.assert_allclose(rotation, np.diag([-1, -1, 1]))
    registration = estimate_ct_registration(vertices, [20, 0, -300], rotation, mattress_mm=0,
                                            target_distance_from_head_mm=500, mesh_sha256="x", target_info={}, pose_info=pose)
    np.testing.assert_allclose(registration["anchor_ct_lps_mm"], [20, -110, -300])   # anterior (LPS y < 0)


def test_room_reference_puts_the_target_at_the_isocentre():
    registration = {"anchor_ct_lps_mm": A_C, "anchor_table_mm": A_T, "rotation_table_from_ct": np.eye(3).tolist()}
    events = [event(M0, "a", ka=0.), event(M1, "b", ka=1.), event([300, 35, 185], "c", ka=5.)]
    reference = estimate_room_reference(registration, P_C, events)
    assert reference["estimate"]["reference_event_uid"] == "b"          # first event with kerma
    assert reference["readings_long_lat_height_mm"] == M1
    np.testing.assert_allclose(reference["origin_isocenter_mm"], [0, 35, 500])
    f = frames(origin=reference["origin_isocenter_mm"])
    f = RegisteredFrames(f.table_from_ct, np.asarray(M1, float), f.reference_origin_i, S, {})
    np.testing.assert_allclose(f.table_from_ct[:3, :3] @ P_C + f.table_from_ct[:3, 3] + f.table_origin_i(events[1]), 0)
    assert choose_reference_event(events, "max_kerma")[1].uid == "c"
    assert choose_reference_event(events, "2")[1].uid == "b"
    with pytest.raises(InputError):
        choose_reference_event(events, "1")


def test_target_selection_order():
    vertices = cylinder()
    chest_lps = np.array([20., -10., -300.])
    metadata = {"anatomical_anchors_ras_mm": {"chest": (chest_lps*[-1, -1, 1]).tolist()}}
    target, info = select_target(vertices, metadata, ["Heart"])
    np.testing.assert_allclose(target, chest_lps)
    assert info["target_source"] == "ct_anatomical_chest"
    shoulder = {"anchor_ct_lps_mm": [20, 110, -100], "anchor_table_mm": [0, -39, -400],
                "rotation_table_from_ct": np.eye(3).tolist(), "anchor_definition": "posterior_shoulder_midline"}
    target, info = select_target(vertices, {}, [], shoulder)
    np.testing.assert_allclose(target, [20, 0, -280], atol=1e-9)        # 180 mm caudal, mid-depth
    assert info["target_source"].startswith("cardiac_level_proxy")
    target, info = select_target(vertices, metadata, ["Pelvis"])
    assert info["target_source"] == "mesh_bounds_center_geometric_only"
    np.testing.assert_allclose(mesh_vertices_lps([[1, 2, 3]], "RAS"), [[-1, -2, 3]])


def test_table_calibration_from_single_axis_moves():
    events = [event([420, -10, 155], "0"), event([520, -10, 155], "1"), event([520, 40, 155], "2"),
              event([520, 40, 135], "3")]
    profile = {"motion_matrix_i_per_mm": [[1, 0, 0], [0, 0, -1], [0, None, 0]], "motion_matrix_status": {}}
    moves = [parse_move("2:izquierda:100"), parse_move("3:pies"), parse_move("4:arriba:20")]
    result = calibrated_profile(profile, events, moves)
    np.testing.assert_allclose(np.asarray(result["motion_matrix_i_per_mm"], dtype=float),
                               # columns: long +100 -> left 100; lat +50 -> feet; height -20 -> up 20
                               [[1, 0, 0], [0, 0, 1], [0, -1, 0]])
    assert result["motion_matrix_status"]["lateral (DCM 113752)"] == "measured_direction_unit_gain"
    assert result["motion_matrix_status"]["longitudinal (DCM 113751)"] == "measured_rdsr_moves"
    assert np.isfinite(motion_matrix(result)).all()
    with pytest.raises(InputError, match="un solo eje"):
        calibrated_profile(profile, [events[0], replace(events[3], uid="x")], [parse_move("2:arriba")])
    with pytest.raises(InputError, match="Sentido"):
        parse_move("2:atras")
