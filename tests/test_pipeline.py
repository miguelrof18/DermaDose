"""Checks with analytic expectations and independent changes to RDSR fixtures."""
from copy import deepcopy
from dataclasses import replace
from importlib.resources import files
import json
import math

import numpy as np
import pydicom
import pytest

from skindose_verifier.cli import compatibility_report, main, run
from skindose_verifier.config import Equipment, Factor
from skindose_verifier.convergence import compare_resolutions
from skindose_verifier.domain import InputError
from skindose_verifier.engine import calculate, normalize
from skindose_verifier.geometry import beam_from_event, event_slices, project
from skindose_verifier.patient import sample_patient
from skindose_verifier.rdsr import find_rdsr, read_rdsr
from skindose_verifier.sensitivity import position_sensitivity
from skindose_verifier.synthetic import make_synthetic_rdsr
from skindose_verifier.validation import compare_with_radimetrics, validate_measurements, validate_map


@pytest.fixture
def case(tmp_path):
    path = make_synthetic_rdsr(tmp_path / "fictional.dcm")
    source = Equipment.load(__import__("importlib").resources.files(
        "skindose_verifier").joinpath("data/equipo_sintetico.json"))
    conf = tmp_path / "equipment.json"
    conf.write_text(json.dumps(source.raw), encoding="utf-8")
    return path, conf, source, read_rdsr(path)


def config_variant(case, tmp_path, change):
    data = deepcopy(case[2].raw)
    change(data)
    path = tmp_path / "modified.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return Equipment.load(path)


def test_rdsr_event_vs_study_totals_and_units(case):
    parsed = case[3]
    assert len(parsed.events) == 2
    assert parsed.warnings == ()
    assert parsed.accumulated_ka_r_gy == pytest.approx(2 * .019834710744)
    assert parsed.accumulated_kap_gy_m2 == pytest.approx(2 * .00019834710744)
    assert parsed.events[0].exposure_mas == pytest.approx(20)
    assert parsed.events[0].filtration_mm_cu == .2
    assert parsed.events[0].original["ka_r_gy"][0]["unit"] == "Gy"
    assert parsed.events[0].table_system == "isocenter"


def test_find_rdsr_in_dicom_folder_and_reports_non_rdsr(case, tmp_path):
    folder = tmp_path / "dicom_export"
    folder.mkdir()
    rdsr = make_synthetic_rdsr(folder / "dose_sr.dcm")
    ds = pydicom.dcmread(rdsr, stop_before_pixels=True)
    ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.12.1"
    ds.Modality = "XA"
    ds.save_as(folder / "image_xa.dcm", enforce_file_format=True)

    found = find_rdsr(folder)

    assert found["rdsr_count"] == 1
    assert found["rdsr_candidates"][0]["path"].endswith("dose_sr.dcm")
    assert found["rdsr_candidates"][0]["event_count"] == 2
    assert any(item["modality"] == "XA" for item in found["non_rdsr_summary"])


def test_compatibility_report_checks_required_fields_and_config(case, tmp_path):
    folder = tmp_path / "dicom_export"
    folder.mkdir()
    make_synthetic_rdsr(folder / "dose_sr.dcm")

    report = compatibility_report(folder, case[1], "skin_dose")

    item = report["compatibility_reports"][0]
    assert item["readable_by_parser"] is True
    assert item["configuration_device_match"] is True
    assert item["events_with_missing_required_fields"] == 0
    assert item["events_missing_beam_quality"] == 0
    assert item["compatible_for_stage_without_full_dose_run"] is True


@pytest.mark.parametrize("code,new_value,new_unit,expected", [
    ("113738", "19.834710744", "mGy", .019834710744),
    ("122130", "1.9834710744", "Gy.cm2", .00019834710744),
    ("113750", "110", "cm", 1100),
    ("112011", "-35", "deg", -35),
])
def test_unit_conversions_are_not_text_dependent(case, code, new_value, new_unit, expected):
    ds = pydicom.dcmread(case[0])
    event = ds.ContentSequence[0].ContentSequence[0]
    target = next(n for n in event.ContentSequence if n.ConceptNameCodeSequence[0].CodeValue == code)
    target.ConceptNameCodeSequence[0].CodeMeaning = "Intentionally wrong/free text"
    target.MeasuredValueSequence[0].NumericValue = new_value
    target.MeasuredValueSequence[0].MeasurementUnitsCodeSequence[0].CodeValue = new_unit
    ds.save_as(case[0], enforce_file_format=True)
    parsed = read_rdsr(case[0]).events[0]
    field = {"113738": "ka_r_gy", "122130": "kap_gy_m2", "113750": "sid_mm",
             "112011": "primary_deg"}[code]
    assert getattr(parsed, field) == pytest.approx(expected)


def test_unknown_unit_and_missing_data_fail_closed(case, tmp_path):
    ds = pydicom.dcmread(case[0])
    event = ds.ContentSequence[0].ContentSequence[0]
    target = next(n for n in event.ContentSequence if n.ConceptNameCodeSequence[0].CodeValue == "113738")
    target.MeasuredValueSequence[0].MeasurementUnitsCodeSequence[0].CodeValue = "oz"
    ds.save_as(case[0], enforce_file_format=True)
    with pytest.raises(InputError, match="unidad desconocida"):
        read_rdsr(case[0])
    missing = make_synthetic_rdsr(tmp_path / "missing.dcm", omit="113789")
    with pytest.raises(InputError, match="field_width_mm"):
        calculate(read_rdsr(missing).events, case[2])


def test_ka_absent_never_replaced_by_accumulated_total(case, tmp_path):
    path = make_synthetic_rdsr(tmp_path / "no_ka.dcm", omit="113738")
    report = read_rdsr(path)
    assert report.events[0].ka_r_gy is None
    with pytest.raises(InputError, match="ka_r_gy"):
        calculate(report.events, case[2])


def test_missing_beam_quality_blocks_calibrated_result(case):
    ev = replace(case[3].events[0], filtration_mm_cu=None)
    assert calculate([ev], case[2], "free_air").summary["psd_point_gy"] > 0
    with pytest.raises(InputError, match="kVp y material/espesor"):
        calculate([ev], case[2], "skin_dose")


def test_duplicate_uid_and_wrong_device_rejected(case):
    ds = pydicom.dcmread(case[0])
    event2 = ds.ContentSequence[0].ContentSequence[1]
    event2.ContentSequence[1].UID = ds.ContentSequence[0].ContentSequence[0].ContentSequence[1].UID
    ds.save_as(case[0], enforce_file_format=True)
    with pytest.raises(InputError, match="repetido"):
        read_rdsr(case[0])
    with pytest.raises(InputError, match="no coincide"):
        case[2].require_matching_device("Other", "Other")


def test_biplane_requires_characterized_plane_and_keeps_separate_events(case, tmp_path):
    make_synthetic_rdsr(case[0], biplane=True)
    report = read_rdsr(case[0])
    assert [e.plane for e in report.events] == ["113620", "113621"]
    with pytest.raises(InputError, match="sin geometría"):
        equipment = config_variant(case, tmp_path, lambda d: d["geometry"].update(
            planes={"113620": {}}))
        calculate(report.events, equipment)
    equipment = config_variant(case, tmp_path, lambda d: d["geometry"].update(
        planes={"113620": {}, "113621": {}}))
    assert calculate(report.events, equipment).summary["event_count"] == 2


def test_plane_beam_rectangular_mask_and_analytic_inverse_square(case):
    ev = case[3].events[0]
    equipment = case[2]
    skin = sample_patient(equipment.patient)
    projection = project(ev, equipment, skin)
    # Source at z=-700, skin at z=-100: distance 600 mm.
    assert projection.distance_mm[skin.shape[1] * (skin.shape[0] // 2) + skin.shape[1] // 2] == pytest.approx(600.01, rel=.001)
    # Detector 200 mm field -> skin half-width 200/2 * 600/1100=54.545 mm.
    x = np.abs(skin.points_mm[:, 0])
    y = np.abs(skin.points_mm[:, 1])
    expected_mask = (x <= 100 * 600/1100) & (y <= 100 * 600/1100)
    np.testing.assert_array_equal(projection.visible_mask, expected_mask)
    result = calculate([ev], equipment, "free_air")
    assert result.summary["psd_point_gy"] == pytest.approx(
        ev.ka_r_gy * (550/np.min(projection.distance_mm[expected_mask]))**2)
    assert np.count_nonzero(result.total_gy) == expected_mask.sum()


def test_distance_scale_energy_conservation_and_stage_factors(case):
    ev = case[3].events[0]
    equip = case[2]
    air = calculate([ev], equip, "free_air")
    cal = calculate([ev], equip, "calibrated_air")
    absorbed = calculate([ev], equip, "skin_dose")
    mask = air.total_gy > 0
    np.testing.assert_allclose(cal.total_gy[mask], air.total_gy[mask] * .8)
    np.testing.assert_allclose(absorbed.total_gy[mask], cal.total_gy[mask] * 1.25 * 1.05)
    np.testing.assert_allclose(absorbed.total_gy[~mask], 0)
    assert air.summary["psd_point_gy"] > 0
    assert absorbed.summary["psd_neighbourhood_1cm2"]["dose_gy"] <= absorbed.summary["psd_point_gy"]
    metadata = absorbed.summary["factor_metadata"]["default"]
    assert metadata["bsf"]["provenance"]
    assert metadata["bsf"]["identifier"] == absorbed.summary["calibration_versions"]["bsf"]
    assert metadata["bsf"]["representation"] in ("constant", "tabulated_interpolated")


def test_double_source_skin_distance_gives_approximately_quarter_kerma(case):
    ev = case[3].events[0]
    close = calculate([ev], case[2], "free_air")
    far = calculate([replace(ev, sod_mm=1300, sid_mm=1500)], case[2], "free_air")
    center = int(np.argmax(close.total_gy))
    assert far.total_gy[center] / close.total_gy[center] == pytest.approx(.25, rel=1e-4)


def test_table_transmission_only_when_ray_intersects_table(case, tmp_path):
    equipment = config_variant(case, tmp_path, lambda d: d["geometry"].update(table_plane_z_mm=-50))
    ev = case[3].events[0]
    air = calculate([ev], equipment, "free_air")
    calibrated = calculate([ev], equipment, "calibrated_air")
    np.testing.assert_allclose(air.total_gy, calibrated.total_gy)


def test_rotations_have_unit_axes_and_move_source(case):
    ev = case[3].events[0]
    axis = beam_from_event(ev, case[2])
    turn = beam_from_event(replace(ev, primary_deg=90), case[2])
    assert axis.source_mm == pytest.approx([0, 0, -700])
    assert turn.source_mm == pytest.approx([-700, 0, 0])
    assert turn.direction == pytest.approx([1, 0, 0])
    assert np.dot(turn.direction, turn.transverse_x) == pytest.approx(0)
    other = beam_from_event(replace(ev, secondary_deg=90), case[2])
    assert other.source_mm == pytest.approx([0, 700, 0])


def test_table_translation_moves_footprint_in_patient_coordinates(case):
    ev = case[3].events[0]
    reference = calculate([ev], case[2])
    moved = calculate([replace(ev, table_mm=(40, 0, 0))], case[2])
    peak1 = reference.summary["psd_point_patient_mm"][0]
    peak2 = moved.summary["psd_point_patient_mm"][0]
    assert peak2 < peak1
    assert not np.array_equal(reference.total_gy, moved.total_gy)


def test_no_skin_intersection_is_explicitly_reported(case):
    ev = case[3].events[0]
    with pytest.raises(InputError, match="Ningún rayo"):
        calculate([replace(ev, table_mm=(3000, 0, 0))], case[2])
    result = calculate([ev, replace(ev, uid="fictional-other", table_mm=(3000, 0, 0))], case[2])
    assert result.summary["events_outside_model"] == ["fictional-other"]
    assert result.summary["events"][1]["maximum_patient_point_mm"] is None


def test_ellipse_sample_and_front_face(case, tmp_path):
    equipment = config_variant(case, tmp_path, lambda d: d.update(patient={
        "type": "ellipse", "semi_axis_x_mm": 170, "semi_axis_z_mm": 110,
        "length_mm": 500, "origin_on_table_mm": [0, 0, 0], "resolution_mm": 5}))
    skin = sample_patient(equipment.patient)
    assert np.sum(skin.area_mm2) == pytest.approx(500 * 894.7, rel=.01)
    mask = project(case[3].events[0], equipment, skin).visible_mask
    assert np.any(mask)
    assert np.all(skin.points_mm[mask, 2] < 0)
    assert calculate(case[3].events, equipment, "skin_dose").summary["psd_point_gy"] > 0


def test_event_order_does_not_affect_sum_and_map_is_in_patient_frame(case):
    evs = case[3].events
    a = calculate(evs, case[2], "free_air")
    b = calculate(list(reversed(evs)), case[2], "free_air")
    np.testing.assert_array_equal(a.total_gy, b.total_gy)
    np.testing.assert_allclose(a.total_gy, np.sum(a.event_maps_gy, axis=0))
    assert a.summary["psd_point_gy"] > max(x["maximum_gy"] for x in a.summary["events"])
    contributions = a.summary["global_psd_event_contributions"]
    assert [item["rank"] for item in contributions] == [1, 2]
    assert sum(item["dose_at_global_psd_point_gy"] for item in contributions) == pytest.approx(a.summary["psd_point_gy"])


def test_sweep_is_blocked_without_temporal_profile(case, tmp_path):
    ev = replace(case[3].events[0], primary_end_deg=20)
    with pytest.raises(InputError, match="perfil temporal"):
        calculate([ev], case[2])
    equipment = config_variant(case, tmp_path, lambda d: d.update(allow_uniform_sweep=True,
                                                                   sweep_samples=10))
    pieces = event_slices(ev, equipment)
    assert len(pieces) == 10
    assert sum(p.ka_r_gy for p in pieces) == pytest.approx(ev.ka_r_gy)
    assert calculate([ev], equipment).summary["events"][0]["uniform_sweep_assumed"]


def test_factor_interpolation_and_no_extrapolation(case, tmp_path):
    f = Factor.from_dict("cf", {"provenance": "phantom measurement", "identifier": "curve-v1",
                                 "axes": {"kvp": [60, 100], "filtration_mm_cu": [0, 1]},
                                 "grid": [[1, 2], [3, 4]]})
    assert f.evaluate(kvp=80, filtration_mm_cu=.5) == pytest.approx(2.5)
    with pytest.raises(InputError, match="rango calibrado"):
        f.evaluate(kvp=101, filtration_mm_cu=.5)
    with pytest.raises(InputError, match="doble conteo"):
        config_variant(case, tmp_path, lambda d: d["factors"].update(table_total=d["factors"]["table_normal"]))


def test_mode_specific_calibration_is_selected_and_missing_mode_rejected(case, tmp_path):
    ev = case[3].events[0]
    base = calculate([ev], case[2], "calibrated_air").total_gy
    equipment = config_variant(case, tmp_path, lambda d: d.update(factors_by_mode={
        "44491008": {"cf": {"value": 1.5, "provenance": "fictional", "identifier": "FLUORO-CF"}}}))
    np.testing.assert_allclose(calculate([ev], equipment, "calibrated_air").total_gy,
                               base * 1.5)
    with pytest.raises(InputError, match="sin calibración"):
        calculate([replace(ev, mode="different_acquisition_mode")], equipment, "calibrated_air")


def test_table_angle_uses_each_ray(case):
    geo = project(case[3].events[0], case[2], sample_patient(case[2].patient))
    assert geo.table_angle_deg.min() > 0
    assert geo.table_angle_deg.max() > geo.table_angle_deg.min()


def test_refining_grid_converges_on_synthetic_case(case):
    comparison = compare_resolutions(case[3].events, case[2], "skin_dose", [10, 5, 2.5])
    assert comparison["clinical_validation"] is False
    assert abs(comparison["relative_point_psd_change_last_refinement"]) < .003
    assert comparison["point_displacement_mm_last_refinement"] < 10


def test_independent_yield_vs_derived_kap_consistency(case):
    check = calculate([case[3].events[0]], case[2]).summary["events"][0]["irp_verification"]
    assert check["tube_output"]["ka_r_gy"] == pytest.approx(.0003 * 20 * (1000/550)**2)
    assert check["tube_output"]["relative_difference_to_reported"] == pytest.approx(0, abs=1e-8)
    assert check["kap_area"]["ka_axis_gy"] == pytest.approx(.019834710744, rel=1e-8)
    assert "consistency_only" in check["kap_area"]["classification"]


def test_real_measurements_required_and_reference_medium_enforced(case, tmp_path):
    summary = calculate(case[3].events, case[2], "skin_dose").summary
    with pytest.raises(InputError, match="medio dosimétrico"):
        compare_with_radimetrics(summary, .04, .03, "water")
    csv = tmp_path / "measurements.csv"
    csv.write_text("test_id,rdsr_path,config_path,measured_max_gy,standard_uncertainty_gy\n"
                   f"case1,{case[0].name},{case[1].name},0.035,0.002\n", encoding="utf-8")
    output = validate_measurements(csv, "skin_dose")
    assert output["validated"] is False
    assert output["cases"][0]["predicted_gy"] == pytest.approx(summary["psd_point_gy"])


def test_cli_outputs_auditable_and_machine_readable(case, tmp_path, capsys):
    output = run(case[0], case[1], "skin_dose", tmp_path/"output")
    assert output["psd_point_gy"] == pytest.approx(.035258043529, rel=1e-7)
    paths = output["files"]
    for p in paths.values():
        assert __import__("pathlib").Path(p).exists()
    summary = json.loads(__import__("pathlib").Path(paths["summary"]).read_text(encoding="utf-8"))
    assert summary["audit"]["rdsr_sha256"] == case[3].sha256
    assert summary["status"] == "research_unvalidated_no_clinical_use"
    with np.load(paths["numeric_map"]) as saved:
        assert saved["dose_gy"].shape == (72, 72)
        assert saved["event_dose_gy"].shape == (2, 72, 72)
        np.testing.assert_array_equal(saved["event_dose_gy"].sum(axis=0), saved["dose_gy"])
        area = summary["surface_dose_area_approx"]
        assert area["irradiated_area_cm2"] >= area["area_at_least_10_percent_peak_cm2"]
        assert area["area_at_least_10_percent_peak_cm2"] >= area["area_at_least_50_percent_peak_cm2"]
        assert area["area_at_least_50_percent_peak_cm2"] >= area["area_at_least_90_percent_peak_cm2"] > 0
    assert __import__("pathlib").Path(paths["interactive_3d"]).stat().st_size > 1_000_000
    assert main(["inspect", str(case[0])]) == 0
    assert "PatientName" not in capsys.readouterr().out


def _mesh_equipment(case, tmp_path):
    from shutil import copyfile
    asset = tmp_path / "adult_female.stl"
    copyfile(files("skindose_verifier").joinpath("data/pyskindose_adult_female.stl"), asset)
    equipment = config_variant(case, tmp_path, lambda d: d.update(patient={
        "type": "mesh", "mesh_path": str(asset), "mesh_unit": "cm",
        "mesh_axes_matrix": [[1, 0, 0], [0, 0, 1], [0, -1, 0]],
        "mesh_origin_mm": [0, 350, -100], "thickness_mm": 220,
        "origin_on_table_mm": [0, 0, 0], "resolution_mm": 5}))
    return equipment


def test_community_anthropomorphic_mesh_and_units(case, tmp_path):
    equipment = _mesh_equipment(case, tmp_path)
    surface = sample_patient(equipment.patient)
    assert surface.kind == "mesh"
    assert surface.mesh.is_watertight and surface.mesh.is_volume
    assert surface.shape == (26756,)
    assert surface.mesh.bounds[0] == pytest.approx([-777.8106, -1256.2302, -100], rel=1e-5)
    assert surface.mesh.bounds[1] == pytest.approx([777.8106, 350, 179.4599], rel=1e-5)
    assert np.all(surface.area_mm2 > 0)


def test_mesh_dose_hits_first_surface_and_accumulates(case, tmp_path):
    equipment = _mesh_equipment(case, tmp_path)
    surface = sample_patient(equipment.patient)
    ev = case[3].events[0]
    projection = project(ev, equipment, surface)
    irradiated = np.flatnonzero(projection.visible_mask)
    assert 50 < len(irradiated) < 1000
    assert np.all(surface.normals[irradiated, 2] < 0)
    source = projection.beam.source_mm
    first = surface.mesh.ray.intersects_first(np.broadcast_to(source, (len(irradiated), 3)),
                                              surface.points_mm[irradiated]-source)
    np.testing.assert_array_equal(first, irradiated)
    output = calculate(case[3].events, equipment, "skin_dose")
    assert output.summary["psd_point_gy"] == pytest.approx(.03353185114, rel=1e-6)
    assert 0 < output.summary["psd_neighbourhood_1cm2"]["dose_gy"] <= output.summary["psd_point_gy"]
    assert output.summary["mesh_geometry"]["face_count"] == 26756
    np.testing.assert_allclose(output.total_gy, sum(output.event_maps_gy))


def test_mesh_own_shadow_blocks_hidden_front_facing_triangles(case, tmp_path):
    equipment = _mesh_equipment(case, tmp_path)
    skin = sample_patient(equipment.patient)
    ev = replace(case[3].events[0], primary_deg=45, field_width_mm=500,
                 field_height_mm=600)
    hit = project(ev, equipment, skin)
    beam = hit.beam
    delta = skin.points_mm-beam.source_mm
    depth = delta@beam.direction
    front = np.einsum("ij,ij->i", skin.normals, -delta)>1e-9
    geometrically_inside = (depth>0)&(depth<=ev.sid_mm)&front
    geometrically_inside &= np.abs(delta@beam.transverse_x) <= ev.field_width_mm*depth/(2*beam.reference_distance_mm)
    geometrically_inside &= np.abs(delta@beam.transverse_y) <= ev.field_height_mm*depth/(2*beam.reference_distance_mm)
    # The arm/torso silhouette has faces oriented towards the source that
    # still lie behind another face. First-hit geometry must mask those.
    assert geometrically_inside.sum()-hit.visible_mask.sum() > 100
    assert hit.visible_mask.sum() > 100


def test_mesh_dose_rejects_missing_units_and_open_skin(case, tmp_path):
    equipment = _mesh_equipment(case, tmp_path)
    wrong = deepcopy(equipment.patient)
    wrong.pop("mesh_unit")
    with pytest.raises(InputError, match="STL no guarda unidades"):
        sample_patient(wrong)
    import trimesh
    mesh = trimesh.creation.box(extents=(80, 80, 80))
    mesh.update_faces(np.arange(len(mesh.faces)-1))
    path = tmp_path / "open_skin.stl"
    mesh.export(path)
    wrong = {**equipment.patient, "mesh_path": str(path), "mesh_unit": "mm"}
    with pytest.raises(InputError, match="Malla abierta"):
        sample_patient(wrong)


def test_mesh_convergence_requires_new_mesh_not_a_new_metadata_value(case, tmp_path):
    equipment = _mesh_equipment(case, tmp_path)
    surface = sample_patient(equipment.patient)
    changed = sample_patient({**equipment.patient, "resolution_mm": 1})
    np.testing.assert_array_equal(surface.points_mm, changed.points_mm)


def test_patient_position_scenarios_recalculate_without_mutating_equipment(case):
    original = deepcopy(case[2].raw)
    report = position_sensitivity(case[3].events, case[2], "skin_dose", 5.)
    assert report["base_psd_point_gy"] == pytest.approx(.035258043529, rel=1e-7)
    assert len(report["scenarios"]) == 6
    assert {tuple(row["offset_mm"]) for row in report["scenarios"]} == {
        (-5., 0., 0.), (5., 0., 0.), (0., -5., 0.), (0., 5., 0.), (0., 0., -5.), (0., 0., 5.)}
    assert all(row["status"] == "calculated" for row in report["scenarios"])
    assert report["psd_range_valid_scenarios_gy"][0] < report["psd_range_valid_scenarios_gy"][1]
    assert case[2].raw == original
    assert report["clinical_validation"] is False
    with pytest.raises(InputError, match="desplazamiento"):
        position_sensitivity(case[3].events, case[2], "skin_dose", 0)


def test_mesh_demo_packages_3d_faces_and_mesh_hash(case, tmp_path):
    assert main(["demo", "--model", "mesh", "--out", str(tmp_path/"demo"),
                 "--stage", "skin_dose"]) == 0
    directory = tmp_path / "demo"
    result = json.loads((directory/"resultados.json").read_text(encoding="utf-8"))
    assert result["audit"]["mesh_sha256"] == "2edd0de03344a1081f2a6567bbcca87b965fbbadd822756e05830cd570cdcb3b"
    assert result["psd_point_gy"] == pytest.approx(.034598271523, rel=1e-6)
    assert result["psd_neighbourhood_1cm2"]["sampled_area_mm2"] == pytest.approx(92.68549235, rel=1e-4)
    assert result["mesh_geometry"]["face_count"] == 12148
    with np.load(tmp_path/"demo/mapa_dosis.npz") as contents:
        assert contents["faces"].shape == (12148, 3)
        assert contents["vertices_mm"].shape[1] == 3
        assert contents["dose_gy"].shape == (12148,)
        assert contents["event_dose_gy"].shape == (2, 12148)
        np.testing.assert_array_equal(contents["event_dose_gy"].sum(axis=0), contents["dose_gy"])
    from html.parser import HTMLParser

    class JsonExtractor(HTMLParser):
        reading = False
        payload = ""

        def handle_starttag(self, tag, attrs):
            self.reading = tag == "script" and ("id", "viewer-data") in attrs

        def handle_data(self, data):
            if self.reading:
                self.payload += data

        def handle_endtag(self, tag):
            if tag == "script":
                self.reading = False

    html = (tmp_path/"demo/mapa_3d_interactivo.html").read_text(encoding="utf-8")
    assert "<script src=" not in html and "Planos de corte móviles" in html
    parser = JsonExtractor()
    parser.feed(html)
    viewer = json.loads(parser.payload)
    assert len(viewer["faces"]) == 12148
    assert len(viewer["doses"]) == 3
    assert max(viewer["doses"][0]) == pytest.approx(result["psd_point_gy"], abs=1e-11)
    equipment = Equipment.load(directory/"equipo_sintetico.json")
    surface = sample_patient(equipment.patient)
    projection = project(case[3].events[0], equipment, surface)
    illuminated = np.flatnonzero(projection.visible_mask)
    assert len(illuminated) > 25
    assert np.all(surface.normals[illuminated, 2] < 0)
    first = surface.mesh.ray.intersects_first(
        np.broadcast_to(projection.beam.source_mm, (len(illuminated), 3)),
        surface.points_mm[illuminated]-projection.beam.source_mm)
    np.testing.assert_array_equal(first, illuminated)


def test_map_comparison_uses_registered_samples_and_checks_quantity(case, tmp_path):
    from pathlib import Path
    out = run(case[0], case[1], "skin_dose", tmp_path/"model")["files"]
    summary = json.loads(Path(out["summary"]).read_text(encoding="utf-8"))
    with np.load(out["numeric_map"]) as f:
        pts = f["patient_points_mm"]
        doses = f["dose_gy"].ravel()
    hot = np.flatnonzero(doses > 0)[:3]
    cold = np.flatnonzero(doses == 0)[:3]
    csv_path = tmp_path / "measured_map.csv"
    with csv_path.open("w", encoding="utf-8") as f:
        f.write("x_mm,y_mm,z_mm,measured_gy,standard_uncertainty_gy,quantity\n")
        for idx in list(hot)+list(cold):
            f.write(",".join(map(str, [pts[idx][0]+.1, pts[idx][1], pts[idx][2],
                                       doses[idx], .001, summary["quantity"]]))+"\n")
    comparison = validate_map(out["summary"], csv_path)
    assert comparison["rmse_gy"] == pytest.approx(0)
    assert comparison["isodose_dice_at_sampled_points"] == 1
    assert comparison["validated"] is False
    with pytest.raises(InputError, match="Registro espacial"):
        validate_map(out["summary"], csv_path, max_distance_mm=.001)
