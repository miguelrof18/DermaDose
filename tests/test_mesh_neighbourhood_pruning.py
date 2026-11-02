import heapq
import math
import numpy as np
import pytest
import trimesh
from skindose_verifier.patient import SkinSurface
from skindose_verifier.metrics import neighbourhood_psd


@pytest.mark.parametrize('pattern', ['uniform', 'localized', 'gradient'])
def test_pruned_search_matches_exhaustive_graph_search(pattern):
    mesh = trimesh.creation.icosphere(subdivisions=2, radius=10.)
    surface = SkinSurface(mesh.triangles_center, mesh.face_normals, mesh.area_faces,
                          (len(mesh.faces),), (np.array([]), np.array([])), ('',''), 'mesh', mesh)
    values = np.ones(len(mesh.faces))
    if pattern == 'localized':
        values = np.where(mesh.triangles_center[:,0] > 8, 3., 0.)
    elif pattern == 'gradient':
        values = np.arange(len(mesh.faces), dtype=float)/len(mesh.faces)
    area_mm2 = 50.
    radius = math.sqrt(area_mm2/math.pi)
    graph = surface.face_graph
    best, best_index, best_area = -1., -1, 0.
    for origin in range(len(values)):
        distance = {origin: 0.}
        queue = [(0., origin)]
        selected = []
        while queue:
            length, node = heapq.heappop(queue)
            if length > radius:
                break
            if length > distance[node]:
                continue
            selected.append(node)
            for j in range(graph.indptr[node],graph.indptr[node+1]):
                neighbour = graph.indices[j]
                candidate = length + graph.data[j]
                if candidate <= radius and candidate < distance.get(neighbour, float('inf')):
                    distance[neighbour] = candidate
                    heapq.heappush(queue, (candidate, neighbour))
        weights = surface.area_mm2[selected]
        area = float(weights.sum())
        if not (.85*area_mm2 <= area <= 1.25*area_mm2):
            continue
        mean = float(np.dot(values[selected], weights)/area)
        if mean > best:
            best, best_index, best_area = mean, origin, area
    actual = neighbourhood_psd(surface, values, area_mm2)
    if best_index < 0:
        assert actual is None
    else:
        assert actual['dose_gy'] == best
        assert actual['sampled_area_mm2'] == best_area
        np.testing.assert_array_equal(actual['patient_point_mm'], surface.points_mm[best_index])


def test_mesh_budget_keeps_full_point_dose_and_marks_optional_metric(tmp_path):
    from copy import deepcopy
    from pathlib import Path
    from skindose_verifier.config import Equipment
    from skindose_verifier.engine import calculate
    from skindose_verifier.rdsr import read_rdsr
    from skindose_verifier.synthetic import make_synthetic_rdsr
    mesh = trimesh.creation.icosphere(subdivisions=3, radius=50.)
    surface = SkinSurface(mesh.triangles_center, mesh.face_normals, mesh.area_faces,
                          (len(mesh.faces),), (np.array([]),np.array([])), ('',''), 'mesh', mesh)
    base = Equipment.load(Path(__file__).parents[1]/'skindose_verifier/data/equipo_sintetico.json')
    raw = deepcopy(base.raw)
    raw['max_neighbourhood_faces'] = 1000
    equipment = Equipment.from_dict(raw)
    events = read_rdsr(make_synthetic_rdsr(tmp_path/'source.dcm')).events
    expected = calculate(events,equipment,surface=surface,detailed=False)
    actual = calculate(events,equipment,surface=surface)
    np.testing.assert_array_equal(actual.total_gy,expected.total_gy)
    assert len(actual.total_gy) == len(mesh.faces)
    assert actual.summary['psd_neighbourhood_1cm2'] is None
    assert actual.summary['psd_neighbourhood_1cm2_status'] == 'not_computed_mesh_budget'
    assert actual.summary['psd_point_gy'] == expected.summary['psd_point_gy']
