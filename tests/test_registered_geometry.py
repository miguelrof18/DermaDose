"""Known coordinates and end-to-end checks on a fictitious mesh and RDSR."""
from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pydicom
import pytest
import trimesh

from skindose_verifier.azurion import prepare_events
from skindose_verifier.config import Equipment
from skindose_verifier.domain import InputError
from skindose_verifier.engine import calculate, normalize
from skindose_verifier.example_runner import main as launch
from skindose_verifier.geometry import beam_from_event, beam_view, intersects_tabletop, project, patient_pose, table_shift
from skindose_verifier.patient import sample_patient
from skindose_verifier.rdsr import nodes, read_rdsr
from skindose_verifier.reference_frames import ct_to_table, marker_in_isocenter, motion_matrix, ROOM_FROM_I, RAS_TO_LPS
from skindose_verifier.synthetic import item, make_synthetic_rdsr


ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "config/mesa_azurion_hemodinamica_1.json"
EQUIPMENT = ROOT / "skindose_verifier/data/equipo_sintetico.json"
# Numerical example in the reference-frame document, not measurements of the Azurion.
M0 = [420, -10, 155]
M1 = [320, 35, 185]
S = [[0, 1, 0], [0, 0, -1], [-1, 0, 0]]


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path, data):
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def case(folder, *, axes="RAS"):
    table = read(TABLE)
    table.update(id="synthetic-marker-test", motion_matrix_i_per_mm=S,
                 reference={"readings_long_lat_height_mm": M0, "marker_table_mm": [95, -5, -225],
                            "marker_isocenter_mm": [95, 120, 275]})
    mesh = trimesh.creation.box(extents=[36, 20, 50])
    mesh.apply_translation([20, 100, -360])
    if axes == "RAS":
        mesh.apply_transform(np.diag([-1, -1, 1, 1]))
    path = folder / "piel.stl"
    mesh.export(path)
    digest = sha256(path.read_bytes()).hexdigest()
    registration = {"anchor_ct_lps_mm": [20, 110, -140], "anchor_table_mm": [0, -35, -280],
                    "rotation_table_from_ct": np.eye(3).tolist(), "mesh_sha256": digest}
    write(path.with_suffix(".json"), {"mesh_file": path.name, "mesh_sha256": digest,
          "coordinate_system": axes, "length_unit": "mm", "frame_of_reference_uid": "2.25.101"})
    data = read(EQUIPMENT)
    data["geometry"].update(table_coordinates="legacy", table_calibration=table)
    data["patient"].update(type="mesh", mesh_path=str(path), ct_registration=registration)
    dicom = make_synthetic_rdsr(folder / "serie.dcm")
    ds = pydicom.dcmread(dicom)
    for event, values in zip(nodes(ds, "113706"), [M0, M1]):
        for code, value in zip(("113751", "113752", "113753"), values):
            event.ContentSequence.append(item("NUM", code, "Table reading", value=value / 10, unit="cm"))
        # Keep the displaced box inside the second exposure for dose accumulation.
        if values == M1:
            for code in ("113788", "113789"):
                nodes(event, code)[0].MeasuredValueSequence[0].NumericValue = "560"
    ds.save_as(dicom, enforce_file_format=True)
    return data, dicom


@pytest.mark.parametrize("axes", ["RAS", "LPS"])
def test_matrices_ray_and_dose_in_ct_coordinates(tmp_path, axes):
    data, path = case(tmp_path, axes=axes)
    # A registered case needs no zero, plane height or orientation from the older model.
    for key in ("table_matrix", "table_origin_mm", "table_plane_z_mm"):
        del data["geometry"][key]
    data["azurion_profile"] = {"name": "philips_azurion", "allow_research_assumptions": True,
                                "reference_table_from_first_event": True}
    eq = Equipment.from_dict(data)
    raw = read_rdsr(path).events
    events = [normalize(ev, eq) for ev in prepare_events(raw, eq)]
    skin = sample_patient(eq.patient)
    ct_point = [20, 110, -360, 1]
    for ev, readings, expected_i, origin in zip(events, [M0, M1], [[0, 90, 0], [45, 60, 100]],
                                                [[0, 125, 500], [45, 95, 600]]):
        assert ev.quality["table_reference"] == "fixed_tabletop_marker_reference"
        view = beam_view(ev, eq)
        assert view["readings_long_lat_height_mm"] == readings
        ct_i = np.asarray(view["M_I_from_C"])
        np.testing.assert_allclose(ct_i @ ct_point, [*expected_i, 1])
        np.testing.assert_allclose(np.asarray(view["M_I_from_T"])[:3, 3], origin)
        np.testing.assert_allclose(ct_i, np.asarray(view["M_I_from_T"]) @ view["M_T_from_C"])
        np.testing.assert_allclose(np.asarray(view["M_C_from_I"]) @ ct_i, np.eye(4))
        np.testing.assert_allclose(ct_i @ [*view["isocenter_patient_mm"], 1], [0, 0, 0, 1])
        focus, irp = map(np.asarray, (view["source_patient_mm"], view["irp_ct_lps_mm"]))
        np.testing.assert_allclose(irp - focus, ev.srp_mm * np.array(view["direction_ct_lps"]))
        np.testing.assert_allclose(np.linalg.norm(np.array(view["detector_patient_mm"]) - focus), ev.sid_mm)
    # Independent frontal solution: focus = (20, 720, -360) in CT. Only the posterior
    # face (y=110) is reached. All its rays cross the board and have the same CF/BSF.
    frontal = calculate([raw[0]], eq, "skin_dose", surface=skin)
    expected = np.zeros(len(skin.points_mm))
    posterior = np.isclose(skin.points_mm[:, 1], 110)
    assert posterior.sum() == 2
    distances = np.linalg.norm(skin.points_mm[posterior] - [20, 720, -360], axis=1)
    expected[posterior] = raw[0].ka_r_gy * (550 / distances)**2 * .8 * 1.25 * 1.05
    np.testing.assert_allclose(frontal.total_gy, expected, rtol=1e-12)
    assert frontal.summary["psd_point_gy"] == pytest.approx(expected.max(), rel=1e-12)
    np.testing.assert_array_equal(project(events[0], eq, skin).table_intersection_mask, posterior)
    oblique = beam_from_event(replace(events[0], primary_deg=10, secondary_deg=15), eq)
    assert oblique.detector_mm[0] > 0 and oblique.detector_mm[1] > 0


def test_rotated_patient_and_oblique_moving_table_against_convex_box(tmp_path):
    data, path = case(tmp_path)
    a = np.deg2rad(28)
    rotation = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    data["patient"]["ct_registration"]["rotation_table_from_ct"] = rotation.tolist()
    # An arbitrary room origin must cancel in all coordinates relative to isocentre.
    data["geometry"]["isocenter_mm"] = [13, -25, 47]
    eq = Equipment.from_dict(data)
    skin = sample_patient(eq.patient)
    events = read_rdsr(path).events
    result = calculate(events, eq, "free_air", surface=skin)
    normals_i = skin.normals @ rotation.T
    p_t = (skin.points_mm - [20, 110, -140]) @ rotation.T + [0, -35, -280]
    expected_total = np.zeros(len(p_t))
    for event, origin, per_event in zip(events, [[0, 125, 500], [45, 95, 600]], result.event_maps_gy):
        # Closed-form ray angles in I for these primary-only rotations.
        theta = np.deg2rad(event.primary_deg)
        along = np.array([np.sin(theta), -np.cos(theta), 0])
        across = np.array([np.cos(theta), np.sin(theta), 0])
        rays = p_t + origin + 700 * along
        depth = rays @ along
        visible = (normals_i * rays).sum(axis=1) < 0  # convex box
        visible &= np.abs(rays @ across) <= event.field_width_mm * depth / 2200
        visible &= np.abs(rays[:, 2]) <= event.field_height_mm * depth / 2200
        assert visible.sum() >= 2
        dose = np.zeros_like(expected_total)
        dose[visible] = event.ka_r_gy * 550**2 / (rays[visible]**2).sum(axis=1)
        np.testing.assert_allclose(per_event, dose, rtol=1e-12)
        expected_total += dose
    np.testing.assert_allclose(result.total_gy, expected_total, rtol=1e-12)


def test_measured_board_and_missing_references(tmp_path):
    data, path = case(tmp_path)
    data["geometry"]["table_calibration"] = str(TABLE)
    registration = data["patient"].pop("ct_registration")
    eq = Equipment.from_dict(data)
    assert eq.geometry["table_half_length_mm"] == 1580 and eq.geometry["table_half_width_mm"] == 250
    assert eq.geometry["table_calibration"]["thickness_mm"] == 30
    assert eq.geometry["table_calibration"]["radimetrics"]["offsets_mm"] == {
        "lateral": 508, "longitudinal": -70.5, "height": -732}
    data["patient"]["ct_registration"] = registration
    provisional = Equipment.from_dict(data)
    assert provisional.registered_frames is None
    assert provisional.patient["ct_table_registration"]["ct_to_table_status"] == "explicit"
    table = read(TABLE)
    table["reference"] = {"readings_long_lat_height_mm": M0, "marker_table_mm": [0, 0, 0],
                          "marker_isocenter_mm": [0, 0, 0]}
    data["geometry"]["table_calibration"] = table
    with pytest.raises(InputError, match="motion_matrix"):
        Equipment.from_dict(data)
    table["motion_matrix_i_per_mm"] = S
    eq = Equipment.from_dict(data)
    only_iso = read_rdsr(make_synthetic_rdsr(tmp_path / "solo_iso.dcm")).events[0]
    with pytest.raises(InputError, match="113751.113753"):
        normalize(only_iso, eq)
    ev = normalize(read_rdsr(path).events[0], eq)
    for code, value in (("113754", 10), ("113759", M0[0] + 1)):
        with pytest.raises(InputError, match=code):
            eq.validate_event(replace(ev, original={**ev.original, "table_motion": {code: [value]}}), "free_air")
    with pytest.raises(InputError, match="Plano"):
        eq.geometry_for("113621")
    registration["frame_of_reference_uid"] = "2.25.999"
    with pytest.raises(InputError, match="UID"):
        Equipment.from_dict(data)


def test_fiducials_scale_and_marker_reprojection():
    ct = np.array([[1, 1, 2], [8, 1, 2], [2, 3, 4], [1, 2, 3]])
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    offset = [12, -5, 300]
    t = ct @ rotation.T + offset
    fitted, audit = ct_to_table({"fiducials_ct_lps_mm": ct.tolist(), "fiducials_table_mm": t.tolist()})
    np.testing.assert_allclose(fitted[:3, :3], rotation, atol=1e-12)
    np.testing.assert_allclose(fitted[:3, 3], offset)
    assert audit["fit_rms_mm"] < 1e-12
    line = [[0, 0, 0], [1, 0, 0], [4, 0, 0]]
    with pytest.raises(InputError, match="triángulo"):
        ct_to_table({"fiducials_ct_lps_mm": line, "fiducials_table_mm": line})
    target = np.asarray(S, dtype=float) @ np.diag([1.02, .98, 1.])
    trials = [{"delta_long_lat_height_mm": v.tolist(), "delta_isocenter_mm": (target @ v).tolist()}
              for v in np.eye(3)[[2, 0, 1]] * 20]
    np.testing.assert_allclose(motion_matrix({"axis_trials": trials}), target)
    with pytest.raises(InputError, match="eje"):
        motion_matrix({"axis_trials": [trials[0]] * 3})
    q = [95, 120, 275]
    views = []
    for source in ([0, 0, -500], [200, 0, -500]):
        p = np.diag([800, 800, 1]) @ np.column_stack([np.eye(3), -np.array(source)])
        uvw = p @ [*q, 1]
        views.append({"projection_px_from_i_mm": p.tolist(), "pixel_uv": (uvw[:2] / uvw[2]).tolist()})
    position, stats = marker_in_isocenter({"calibrated_views": views})
    np.testing.assert_allclose(position, q, atol=1e-9)
    assert stats["marker_reprojection_rms_px"] < 1e-9
    with pytest.raises(InputError, match="separación"):
        marker_in_isocenter({"calibrated_views": [views[0]] * 2})


def test_bounded_table_interior_grazing_and_head_edge():
    source = np.array([0., 0., -20.])
    delta = np.array([[0, 0, 40], [0, 25, 40], [0, -25, 40], [60, 0, 40], [0, 0, 5]])
    # Board spans y=-100..0 and z=-10..0 in these local room coordinates.
    hits = intersects_tabletop(source, delta, [-10, -100, -10], [10, 0, 0])
    np.testing.assert_array_equal(hits, [True, False, True, False, False])


def test_launcher_exports_registered_lps_map(tmp_path):
    data, _ = case(tmp_path)
    write(tmp_path / TABLE.name, data["geometry"].pop("table_calibration"))
    write(tmp_path / "registro_tc.json", data["patient"].pop("ct_registration"))
    write(tmp_path / "equipo.json", data)
    assert launch([str(tmp_path)]) == 0
    results = tmp_path / "resultados"
    audit = read(results / "seleccion_entrada.json")
    report = read(results / "resultados.json")
    assert audit["registration_status"] == "explicit"
    assert report["table_placement"] == "marker_reference" and report["event_count"] == 2
    assert all(event["irradiated_cells"] > 0 for event in report["events"])
    with np.load(results / "mapa_dosis.npz") as surface:
        assert str(surface["patient_coordinate_system"]) == "CT_LPS_mm"
        np.testing.assert_allclose(surface["dose_gy"], surface["event_dose_gy"].sum(axis=0))
        i = surface["dose_gy"].argmax()
        np.testing.assert_allclose(surface["patient_points_mm"][i], report["psd_point_patient_mm"])
        assert surface["dose_gy"][i] > 0
    assert Equipment.load(results / "configuracion_efectiva.json").registered_frames is not None
    assert (results / "mapa_3d_interactivo.html").is_file()


@pytest.mark.parametrize("axes", ["RAS", "LPS"])
def test_provisional_room_preserves_ct_table_anchor(tmp_path, axes):
    data, path = case(tmp_path, axes=axes)
    data["geometry"]["table_calibration"]["reference"] = {
        "readings_long_lat_height_mm": None, "marker_table_mm": None, "marker_isocenter_mm": None}
    registration = data["patient"]["ct_registration"]
    eq = Equipment.from_dict(data)
    assert eq.registered_frames is None
    patient = eq.patient
    assert patient["ct_registration"] == registration
    evidence = patient["ct_table_registration"]
    assert evidence["table_to_isocenter_status"] == "provisional_legacy_geometry"
    ct = np.asarray(registration["anchor_ct_lps_mm"])
    native = RAS_TO_LPS @ ct if axes == "RAS" else ct
    surface_point = np.asarray(patient["mesh_axes_matrix"]) @ native + patient["mesh_origin_mm"]
    expected_t = np.asarray(registration["anchor_table_mm"])
    for event in read_rdsr(path).events:
        ev = normalize(event, eq)
        shift = table_shift(ev, eq)
        rotation, translation = patient_pose(eq, shift)
        head_edge = shift + [0, eq.geometry["table_half_length_mm"], eq.geometry["table_plane_z_mm"]]
        np.testing.assert_allclose(rotation @ surface_point + translation,
                                   head_edge + ROOM_FROM_I @ expected_t)
    # Reloading the resolved config must not apply the registration twice.
    reloaded = Equipment.from_dict(eq.raw)
    np.testing.assert_allclose(reloaded.patient["mesh_axes_matrix"], patient["mesh_axes_matrix"])
    np.testing.assert_allclose(reloaded.patient["mesh_origin_mm"], patient["mesh_origin_mm"])
    data["patient"]["ct_registration"]["mesh_sha256"] = "invalid"
    with pytest.raises(InputError, match="otra malla"):
        Equipment.from_dict(data)


def test_launcher_preserves_registration_without_room_reference(tmp_path):
    data, _ = case(tmp_path)
    data["geometry"]["table_calibration"]["reference"] = {}
    write(tmp_path / "registro_tc.json", data["patient"].pop("ct_registration"))
    write(tmp_path / "equipo.json", data)
    assert launch([str(tmp_path), "--prepare-only"]) == 0
    audit = read(tmp_path / "resultados" / "seleccion_entrada.json")
    assert audit["registration_status"] == "ct_table_explicit_room_provisional"
    effective = Equipment.load(tmp_path / "resultados" / "configuracion_efectiva.json")
    assert effective.patient["ct_table_registration"]["ct_to_table_status"] == "explicit"
    assert effective.registered_frames is None


def test_legacy_ambiguous_frame_requires_confirmation_and_known_mixed_rejected(tmp_path):
    data, _ = case(tmp_path)
    data["geometry"]["table_calibration"]["reference"] = {}
    patient = data["patient"]
    meta_path = Path(patient["mesh_path"]).with_suffix(".json")
    meta = read(meta_path)
    ct = {"mixed_frames": True}
    meta["patient_position_evidence"] = {"ct": ct}
    write(meta_path, meta)
    with pytest.raises(InputError, match="Metadatos antiguos"):
        Equipment.from_dict(data)
    patient["confirmed_legacy_frame_metadata_mesh_sha256"] = patient["ct_registration"]["mesh_sha256"]
    eq = Equipment.from_dict(data)
    assert eq.patient["ct_table_registration"]["legacy_frame_metadata_operator_confirmed"]
    assert Equipment.from_dict(eq.raw).patient["ct_table_registration"]["legacy_frame_metadata_operator_confirmed"]
    ct["frame_of_reference_uids"] = ["2.25.101", "2.25.102"]
    write(meta_path, meta)
    with pytest.raises(InputError, match="marcos DICOM distintos"):
        Equipment.from_dict(data)
    ct.update(mixed_frames=False, frame_of_reference_uids=["2.25.101"], missing_frame_of_reference_uids=1)
    write(meta_path, meta)
    assert Equipment.from_dict(data).patient["ct_table_registration"]["frame_of_reference_status"] == "incomplete_dicom_uids"


def test_slicer_frame_metadata_distinguishes_missing_from_distinct():
    import importlib.util
    from types import SimpleNamespace
    spec = importlib.util.spec_from_file_location(
        "slicer_patient_position", ROOT / "slicer/DermaDoseSurface/patient_position.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    h = SimpleNamespace(FrameOfReferenceUID="2.25.101")
    missing = SimpleNamespace()
    other = SimpleNamespace(FrameOfReferenceUID="2.25.102")
    ct = module.ct_geometry(np.eye(4), [h, missing])
    assert not ct["mixed_frames"]
    assert ct["frame_of_reference_uids"] == ["2.25.101"]
    assert ct["missing_frame_of_reference_uids"] == 1
    assert module.ct_geometry(np.eye(4), [h, other])["mixed_frames"]
