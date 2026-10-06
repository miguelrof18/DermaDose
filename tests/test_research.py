"""Analytic, topology, invariance and stochastic checks for the research phase."""
from copy import deepcopy
from dataclasses import replace
from importlib.resources import files

import numpy as np
import pytest
import trimesh

from skindose_verifier.config import Equipment, Factor
from skindose_verifier.convergence import compare_resolutions
from skindose_verifier.domain import InputError
from skindose_verifier.engine import calculate
from skindose_verifier.geometry import project, rot_y
from skindose_verifier.mesh_processing import refine_mesh
from skindose_verifier.metrics import surface_metrics, weighted_dose_metrics, dose_area_summary
from skindose_verifier.patient import sample_patient, SkinSurface
from skindose_verifier.rdsr import read_rdsr
from skindose_verifier.synthetic import make_synthetic_rdsr
from skindose_verifier.uncertainty import propagate_parameters, draw_parameters


@pytest.fixture
def setup(tmp_path):
    equipment = Equipment.load(files("skindose_verifier").joinpath("data/equipo_sintetico.json"))
    events = read_rdsr(make_synthetic_rdsr(tmp_path/"synthetic.dcm")).events
    return equipment, events


def mesh_equipment(equipment, mesh, tmp_path):
    path = tmp_path/"surface.ply"; mesh.export(path)
    cfg = deepcopy(equipment.raw)
    cfg["patient"] = {"type": "mesh", "mesh_path": str(path), "mesh_unit": "mm",
                       "mesh_axes_matrix": np.eye(3).tolist(), "mesh_origin_mm": [0, 0, 0],
                       "resolution_mm": 10, "thickness_mm": 200}
    return Equipment.from_dict(cfg)


@pytest.mark.parametrize("roi", [None, [[-40, -40, -40], [0, 40, 40]]])
def test_conforming_refinement_conserves_area_volume_and_parent_faces(roi):
    original = trimesh.creation.icosphere(subdivisions=1, radius=30)
    refined = refine_mesh(original, 6, roi_mm=roi)
    assert refined.is_watertight and refined.is_winding_consistent
    assert refined.area == pytest.approx(original.area, rel=1e-12)
    assert refined.volume == pytest.approx(original.volume, rel=1e-12)
    parent = refined.metadata["parent_face_index"]
    np.testing.assert_allclose(np.bincount(parent, weights=refined.area_faces), original.area_faces, rtol=1e-12)
    assert np.all(refined.area_faces > 0)
    if roi is None:
        assert refined.edges_unique_length.max() <= 6*(1+1e-12)
    with pytest.raises(InputError, match="mesh_max_faces"):
        refine_mesh(original, .1, max_faces=len(original.faces))


def test_off_axis_field_area_is_at_axial_depth(setup):
    eq, events = setup
    s = sample_patient(eq.patient); g = project(events[0], eq, s)
    expected = 200*200*(600/1100)**2/100
    np.testing.assert_allclose(g.field_area_cm2, expected)
    expected_angle = np.degrees(np.arctan2(np.linalg.norm(s.points_mm[:, :2], axis=1), 600))
    np.testing.assert_allclose(g.skin_angle_deg, expected_angle, atol=1e-10)


def test_calibrated_profile_and_angular_bsf_against_closed_formula(setup):
    eq, events = setup; cfg = deepcopy(eq.raw)
    cfg["factors"]["beam_profile"] = {"identifier": "synthetic-profile", "provenance": "analytic fixture",
        "axes": {"field_fraction_x": [-1, 0, 1]}, "grid": [.8, 1, 1.2]}
    cfg["factors"]["bsf"] = {"identifier": "synthetic-bsf", "provenance": "analytic fixture",
        "axes": {"skin_angle_deg": [0, 90]}, "grid": [1, 1.9]}
    corrected = Equipment.from_dict(cfg)
    skin = sample_patient(eq.patient); geo = project(events[0], eq, skin)
    baseline = calculate([events[0]], eq, "calibrated_air")
    actual = calculate([events[0]], corrected, "skin_dose")
    ids = geo.visible_mask
    expected = baseline.total_gy[ids]*(1+.2*geo.field_fraction_xy[ids, 0])*(1+geo.skin_angle_deg[ids]/100)*1.05
    np.testing.assert_allclose(actual.total_gy[ids], expected, rtol=1e-13)


def test_cached_geometry_is_bitwise_equal_and_reduces_ray_work(setup, tmp_path):
    eq, events = setup
    eq = mesh_equipment(eq, trimesh.creation.icosphere(subdivisions=2, radius=100), tmp_path)
    repeated = [replace(events[0], uid=str(i)) for i in range(6)]
    cached = calculate(repeated, eq, cache_geometry=True)
    plain = calculate(repeated, eq, cache_geometry=False)
    np.testing.assert_array_equal(cached.total_gy, plain.total_gy)
    assert cached.summary["performance"]["projection_cache_hits"] == 5
    assert 6*cached.summary["performance"]["candidate_triangle_pairs"] == plain.summary["performance"]["candidate_triangle_pairs"]


def test_sphere_entry_distance_and_no_back_surface(setup, tmp_path):
    eq, events = setup
    errors = []
    for subdivisions in (1, 3):
        eqm = mesh_equipment(eq, trimesh.creation.icosphere(subdivisions=subdivisions, radius=100), tmp_path)
        surface = sample_patient(eqm.patient)
        ev = replace(events[0], field_width_mm=800, field_height_mm=800)
        g = project(ev, eqm, surface); ids = g.visible_mask
        source = g.beam.source_mm
        ray = surface.points_mm[ids]-source
        ray /= np.linalg.norm(ray, axis=1)[:, None]
        dot = ray @ source
        # Independent ray/sphere quadratic, radius exactly 100 mm.
        roots = -dot-np.sqrt(np.maximum(dot**2-(source@source-100**2), 0))
        expected = ev.ka_r_gy*(ev.srp_mm/roots)**2
        actual = calculate([ev], eqm).total_gy[ids]
        errors.append(float(np.mean(np.abs(actual/expected-1))))
        assert np.all(surface.points_mm[ids, 2] < 0)
    assert errors[1] < errors[0]/3


def test_cylinder_closed_side_entry_and_inverse_square(setup, tmp_path):
    eq, events = setup
    mesh = trimesh.creation.cylinder(radius=100, height=300, sections=64)
    mesh.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2, [1, 0, 0]))
    eq = mesh_equipment(eq, mesh, tmp_path)
    ev = events[0]; result = calculate([ev], eq)
    ids = result.total_gy > 0
    points = result.surface.points_mm[ids]
    assert np.all(points[:, 2] < 0)
    expected = ev.ka_r_gy*ev.srp_mm**2/np.sum((points-[0, 0, -700])**2, axis=1)
    np.testing.assert_allclose(result.total_gy[ids], expected, rtol=1e-13)


def test_open_cylinder_end_check_uses_tilted_patient_frame(setup):
    eq, events = setup
    cfg = deepcopy(eq.raw)
    cfg["patient"] = {"type": "ellipse", "semi_axis_x_mm": 100, "semi_axis_z_mm": 100,
                       "length_mm": 400, "resolution_mm": 10, "tilt_x_deg": 60}
    rotated = Equipment.from_dict(cfg)
    with pytest.raises(InputError, match="extremos"):
        project(events[0], rotated, sample_patient(rotated.patient))


def test_sphere_reflection_symmetry_for_opposite_gantry_angles(setup, tmp_path):
    eq, events = setup
    eq = mesh_equipment(eq, trimesh.creation.icosphere(subdivisions=2, radius=100), tmp_path)
    left = calculate([replace(events[0], primary_deg=25)], eq)
    right = calculate([replace(events[0], primary_deg=-25)], eq)
    np.testing.assert_allclose(np.sort(left.total_gy), np.sort(right.total_gy), atol=1e-12)


def test_winding_and_coordinate_frame_invariance(setup, tmp_path):
    eq, events = setup
    mesh = trimesh.creation.icosphere(subdivisions=2, radius=100)
    eq1 = mesh_equipment(eq, mesh, tmp_path)
    result = calculate([events[0]], eq1)
    mesh.invert(); mesh.vertices = np.asarray(mesh.vertices) @ rot_y(90).T
    eq2 = mesh_equipment(eq, mesh, tmp_path)
    eq2.raw["patient"]["mesh_axes_matrix"] = rot_y(-90).tolist()
    second = calculate([events[0]], eq2)
    # File loading may reorder faces: compare distributions and integral.
    np.testing.assert_allclose(np.sort(result.total_gy), np.sort(second.total_gy), rtol=2e-7, atol=1e-12)
    assert second.surface.mesh.volume > 0


def test_area_metrics_do_not_weight_large_and_small_faces_equally():
    dose = np.array([10., 5., 0.]); area = np.array([1., 9., 90.])
    m = weighted_dose_metrics(dose, area)
    assert m["dmean_gy"] == .55
    assert m["d2_percent_gy"] == 5
    assert m["d10_percent_gy"] == 5
    assert weighted_dose_metrics(np.array([]), np.array([]))["dmean_gy"] is None


def test_hotspots_are_connected_and_histogram_is_monotone(setup, tmp_path):
    eq, events = setup
    eq = mesh_equipment(eq, trimesh.creation.icosphere(subdivisions=2, radius=100), tmp_path)
    result = calculate([events[0]], eq)
    m = result.summary["surface_metrics"]
    histogram = m["cumulative_dose_area_histogram"]
    assert np.all(np.diff(histogram["area_at_least_cm2"]) <= 0)
    assert m["hotspots"]["region_count"] >= 1
    assert m["irradiated_surface"]["dmean_gy"] <= m["whole_model_surface"]["dmax_gy"]
    zero = dose_area_summary(result.surface, np.zeros_like(result.total_gy))
    assert zero["area_at_least_90_percent_peak_cm2"] == 0
    bipolar = (np.abs(result.surface.points_mm[:, 2]) > 85).astype(float)
    assert surface_metrics(result.surface, bipolar)["hotspots"]["region_count"] == 2


def test_real_mesh_convergence_retains_geometry_and_improves_sampling(setup, tmp_path):
    eq, events = setup
    eq = mesh_equipment(eq, trimesh.creation.box([180, 180, 200]), tmp_path)
    results = compare_resolutions([events[0]], eq, "free_air", [30, 15, 7.5])
    samples = [v["surface_samples"] for v in results["values"]]
    assert samples[0] < samples[1] < samples[2]
    exact_peak = events[0].ka_r_gy*(550/600)**2
    assert abs(results["values"][-1]["psd_point_gy"]/exact_peak-1) < .0001
    assert results["method"] == "shared_edge_refinement_fixed_anatomy"


def test_parameter_propagation_matches_analytic_uniform_scale_and_seed(setup):
    eq, events = setup
    spec = {"samples": 64, "seed": 123, "provenance": "synthetic test",
            "parameters": {"ka_scale": {"distribution": "uniform", "mean": 1., "standard_deviation": .1}}}
    a = propagate_parameters(events, eq, "free_air", spec)
    b = propagate_parameters(events, eq, "free_air", spec)
    np.testing.assert_array_equal(a.mean_gy, b.mean_gy)
    ratio = np.asarray(a.summary["psd_samples_gy"])/a.summary["base_psd_gy"]
    np.testing.assert_allclose(ratio, np.asarray(a.summary["parameter_samples"])[:, 0], rtol=1e-13)
    assert ratio.mean() == pytest.approx(1, abs=.002)
    assert ratio.std(ddof=1) == pytest.approx(.1, rel=.03)
    assert np.all((0 <= a.probability_exceeding_threshold) & (a.probability_exceeding_threshold <= 1))


def test_zero_uncertainty_and_moving_patient(setup):
    eq, events = setup
    spec = {"samples": 8, "seed": 4, "provenance": "synthetic test", "parameters": {
        "patient_z_mm": {"distribution": "normal", "mean": 5., "standard_deviation": 0.}}}
    result = propagate_parameters(events, eq, "free_air", spec)
    assert result.summary["psd_sample_standard_deviation_gy"] == 0
    assert result.summary["psd_mean_gy"] < result.summary["base_psd_gy"]


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_nonfinite_quality_fails_closed(setup, bad):
    eq, events = setup
    with pytest.raises(InputError):
        calculate([replace(events[0], kvp=bad)], eq, "skin_dose")
    with pytest.raises(InputError):
        calculate([replace(events[0], primary_end_deg=bad)], eq)
