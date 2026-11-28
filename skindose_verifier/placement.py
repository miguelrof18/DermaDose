"""Estimates for the links of the C -> T -> I chain that have not been measured.

The chain itself lives in reference_frames.py (see docs/sistemas_de_referencia.md):

    p_T = R_T<-C p_C + t_T<-C,    t_k = t_0 + S (m_k - m_0),    p_I(k) = p_T + t_k

This module only *produces inputs* for that chain, as plain JSON that is stored
in the effective configuration together with the method and the hypotheses:

* ``estimate_ct_registration``: M_T<-C from the intervention pose of the RDSR,
  a target point of the CT surface, the mattress thickness and the distance
  from the head edge (document, section 02: a_C <-> a_T = (0, -e, -d)).
* ``estimate_room_reference``: t_0 and m_0 assuming that the target point was
  at the isocentre in a reference event (isocentric technique).

Measured data always take precedence; these functions are never called when the
corresponding link is available. NumPy only, so they can be tested in isolation.
"""
from __future__ import annotations

import numpy as np

from .domain import InputError
from .reference_frames import POSE_ROTATIONS, RAS_TO_LPS, ct_to_table, intervention_pose, table_readings

THORACIC_REGIONS = {"thorax", "thoracic region", "heart", "cardiac", "lung", "lungs", "pulmonary", "chest",
                    "tórax", "torax", "corazón", "corazon", "coronary artery", "coronary arteries"}
# Proxy of the cardiac level below the posterior shoulder anchor (XCAT section).
CHEST_BELOW_SHOULDERS_MM = 180.
# Supporting slab around the target level used to find the back on the mattress.
SUPPORT_SLAB_MM = 50.


def mesh_vertices_lps(vertices, coordinate_system: str) -> np.ndarray:
    """STL vertices in CT LPS mm. RAS -> LPS is (x, y, z) -> (-x, -y, z)."""
    vertices = np.asarray(vertices, dtype=float)
    if vertices.ndim != 2 or vertices.shape[1] != 3 or not np.all(np.isfinite(vertices)):
        raise InputError("Vértices de la malla inválidos.")
    if coordinate_system == "RAS":
        return vertices @ RAS_TO_LPS.T
    if coordinate_system == "LPS":
        return vertices.copy()
    raise InputError("La malla debe documentar coordinate_system RAS o LPS.")


def rdsr_readings(event) -> np.ndarray | None:
    """(longitudinal, lateral, height) in mm, whatever table system the parser selected."""
    original = event.original or {}
    keys = ("table_long", "table_lat", "table_height")
    if all(original.get(k) for k in keys):
        return np.array([float(original[k][0]["value"]) for k in keys])
    if event.table_system == "legacy" and event.table_mm is not None:
        return table_readings(event)
    return None


def pose_rotation(events) -> tuple[np.ndarray, dict]:
    """R_T<-C for the coded intervention pose (HFS, FFS, HFP, FFP)."""
    pose = intervention_pose(events)
    if pose["position"] in POSE_ROTATIONS:
        return POSE_ROTATIONS[pose["position"]].copy(), {**pose, "rotation_source": "rdsr_coded_pose"}
    return np.eye(3), {**pose, "position": "HFS",
                       "rotation_source": "assumed_head_first_supine_rdsr_pose_missing"}


def _slab(vertices_t, level_z, slab_mm=SUPPORT_SLAB_MM):
    selected = vertices_t[np.abs(vertices_t[:, 2] - level_z) <= slab_mm]
    if len(selected) < 3:
        raise InputError("La malla no tiene sección suficiente al nivel del punto diana.")
    return selected


def select_target(vertices_lps, metadata: dict, regions=(), registration: dict | None = None):
    """Target point in CT LPS used to place the beam, and its provenance.

    Order: Slicer anchor matching the RDSR target region; Slicer chest anchor
    for thoracic or unknown regions; a cardiac-level proxy below a shoulder
    anchor; otherwise the centre of the mesh bounds (geometric only).
    """
    vertices_lps = np.asarray(vertices_lps, dtype=float)
    regions = {str(r).strip().casefold() for r in regions if r}
    region = next(iter(regions)) if len(regions) == 1 else None
    thoracic = region is None or region in THORACIC_REGIONS
    anchors = metadata.get("anatomical_anchors_ras_mm") or {}
    if not isinstance(anchors, dict):
        anchors = {}
    lower, upper = vertices_lps.min(axis=0), vertices_lps.max(axis=0)

    def inside(point):
        return point.shape == (3,) and np.all(np.isfinite(point)) and np.all(point >= lower) and np.all(point <= upper)

    for key in ([region] if region in anchors else []) + (["chest"] if thoracic and "chest" in anchors else []):
        try:
            candidate = RAS_TO_LPS @ np.asarray(anchors[key], dtype=float)
        except (TypeError, ValueError):
            continue
        if inside(candidate):
            return candidate, {"target_source": f"ct_anatomical_{key}", "target_region": region}
    shoulder_anchor = registration is not None and (
        registration.get("anchor_definition") == "posterior_shoulder_midline"
        # Exports of DermaDose Surface before anchor_definition was recorded.
        or registration.get("orientation_basis") == "operator_confirmed_aligned_head_first_supine")
    if thoracic and shoulder_anchor and registration.get("anchor_ct_lps_mm") is not None:
        transform, _ = ct_to_table(registration)
        rotation, translation = transform[:3, :3], transform[:3, 3]
        anchor_t = rotation @ np.asarray(registration["anchor_ct_lps_mm"], dtype=float) + translation
        level = anchor_t[2] - CHEST_BELOW_SHOULDERS_MM
        surface_t = vertices_lps @ rotation.T + translation
        try:
            section = _slab(surface_t, level, 10.)
        except InputError:
            low, high = surface_t[:, 2].min(), surface_t[:, 2].max()
            raise InputError(
                "El anclaje de hombros del registro TC–tablero no es coherente con la malla: el nivel cardíaco "
                f"estimado ({CHEST_BELOW_SHOULDERS_MM:g} mm caudal al anclaje) queda fuera de ella. Anclaje en el TC "
                f"(LPS): z = {float(registration['anchor_ct_lps_mm'][2]):.0f} mm; la malla ocupa "
                f"{high - low:.0f} mm, de {high - anchor_t[2]:.0f} mm craneal a {anchor_t[2] - low:.0f} mm caudal al "
                "anclaje. Compruebe en DermaDose Surface que el punto está en la línea media posterior de los hombros "
                "(no en el extremo caudal del TC, frecuente en estudios con pies primero) y la distancia "
                "cabecera–hombros del tablero, o retire registro_tc.json para usar la colocación estimada.") from None
        target_t = np.array([anchor_t[0], (section[:, 1].min() + section[:, 1].max()) / 2, level])
        target = rotation.T @ (target_t - translation)
        return target, {"target_source": f"cardiac_level_proxy_{CHEST_BELOW_SHOULDERS_MM:g}mm_caudal_to_shoulders",
                        "target_region": region}
    return (lower + upper) / 2, {"target_source": "mesh_bounds_center_geometric_only", "target_region": region,
                                 "warning": "Centro geométrico de la malla; no es un órgano ni la región irradiada."}


def estimate_ct_registration(vertices_lps, target_lps, rotation, *, mattress_mm: float,
                             target_distance_from_head_mm: float, mesh_sha256: str,
                             target_info: dict, pose_info: dict) -> dict:
    """M_T<-C with the target over the board midline, its back on the mattress.

    a_C is the target projected onto the supporting surface of the patient at the
    target level; a_T = (0, -e, -d). Returned as a ``registro_tc.json`` object.
    """
    rotation = np.asarray(rotation, dtype=float)
    if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-9) or not np.isclose(np.linalg.det(rotation), 1.):
        raise InputError("La orientación TC–tablero debe ser una rotación propia.")
    for name, value in (("espesor de colchoneta", mattress_mm), ("distancia desde la cabecera", target_distance_from_head_mm)):
        if not np.isfinite(value) or value < 0:
            raise InputError(f"{name} inválido.")
    rotated = np.asarray(vertices_lps, dtype=float) @ rotation.T
    target = rotation @ np.asarray(target_lps, dtype=float)
    support = float(_slab(rotated, target[2])[:, 1].max())   # T: y points down, the back is the largest y
    anchor_ct = rotation.T @ np.array([target[0], support, target[2]])
    anchor_table = np.array([0., -float(mattress_mm), -float(target_distance_from_head_mm)])
    return {
        "mesh_sha256": mesh_sha256,
        "anchor_ct_lps_mm": anchor_ct.tolist(),
        "anchor_table_mm": anchor_table.tolist(),
        "rotation_table_from_ct": rotation.tolist(),
        "anchor_definition": "support_surface_below_target",
        "status": "estimated_automatic",
        "estimate": {
            "method": "target_over_board_midline_back_on_mattress",
            **target_info, "target_ct_lps_mm": np.asarray(target_lps, dtype=float).tolist(),
            "mattress_mm": float(mattress_mm), "target_distance_from_head_edge_mm": float(target_distance_from_head_mm),
            "intervention_pose": pose_info,
            "hypotheses": ["Paciente rígido, con la forma del TC durante el procedimiento.",
                           "Punto diana sobre la línea media del tablero.",
                           "Espalda apoyada sobre la colchoneta comprimida al nivel del punto diana.",
                           "La distancia longitudinal a la cabecera solo afecta a qué rayos cruzan el extremo del tablero."]},
    }


def choose_reference_event(events, reference="kerma_weighted"):
    """Event whose table position defines m_0 for an estimated room reference."""
    usable = [(i, e) for i, e in enumerate(events, 1)
              if rdsr_readings(e) is not None and e.ka_r_gy is not None and e.ka_r_gy > 0]
    if not usable:
        raise InputError("Ningún evento tiene las tres lecturas de mesa (DCM 113751–113753) y kerma positivo.")
    if reference == "first":
        return usable[0]
    if reference == "kerma_weighted":
        # Working position: kerma-weighted median of each reading; the event
        # closest to it (ties: larger kerma) defines m_0.
        readings = np.array([rdsr_readings(e) for _, e in usable])
        weights = np.array([e.ka_r_gy for _, e in usable], dtype=float)
        order = np.argsort(readings, axis=0)
        median = np.array([readings[order[np.searchsorted(np.cumsum(weights[order[:, j]]), weights.sum() / 2), j], j]
                           for j in range(3)])
        distance = np.linalg.norm(readings - median, axis=1)
        best = min(range(len(usable)), key=lambda k: (distance[k], -weights[k]))
        return usable[best]
    if reference == "max_kerma":
        return max(usable, key=lambda item: item[1].ka_r_gy)
    try:
        index = int(reference)
    except (TypeError, ValueError) as exc:
        raise InputError("Evento de referencia: 'kerma_weighted', 'first', 'max_kerma' o un número de evento.") from exc
    for i, event in usable:
        if i == index:
            return i, event
    raise InputError(f"El evento {index} no tiene lecturas de mesa completas o kerma positivo.")


def estimate_room_reference(registration: dict, target_lps, events, *, reference="kerma_weighted",
                            target_info: dict | None = None) -> dict:
    """t_0, m_0 such that the target is at the isocentre in the reference event.

    With R_I<-T = Id: p_I = p_T + t_0 = 0  =>  t_0 = -(R_T<-C a + t_T<-C).
    Returned as a ``table_calibration.reference`` object (status estimated_*).
    """
    transform, _ = ct_to_table(registration)
    target_t = transform[:3, :3] @ np.asarray(target_lps, dtype=float) + transform[:3, 3]
    index, event = choose_reference_event(events, reference)
    m0 = rdsr_readings(event)
    t0 = -target_t
    board_top_below_iso = float(t0[1])
    warnings = []
    if not 30. <= board_top_below_iso <= 400.:
        warnings.append(f"Tablero {board_top_below_iso:.0f} mm bajo el isocentro: fuera del intervalo plausible 30–400 mm.")
    return {
        "readings_long_lat_height_mm": m0.tolist(),
        "origin_isocenter_mm": t0.tolist(),
        "status": "estimated_target_at_isocenter",
        "estimate": {
            "method": "target_at_isocenter_in_reference_event",
            "reference_event_index": index, "reference_event_uid": event.uid,
            "reference_policy": str(reference), **(target_info or {}),
            "target_ct_lps_mm": np.asarray(target_lps, dtype=float).tolist(),
            "target_table_mm": target_t.tolist(),
            "board_top_below_isocenter_mm": board_top_below_iso,
            "board_head_edge_from_isocenter_mm": float(t0[2]),
            "warnings": warnings,
            "hypotheses": ["Técnica isocéntrica: el punto diana estaba en el isocentro en el evento de referencia "
                            "(por defecto, la posición de trabajo: mediana de las lecturas ponderada por Ka,r).",
                           "Tablero nivelado y sin giro (R_I<-T = Id).",
                           "Los demás eventos se desplazan con S (m_k - m_0); no es una calibración de sala."]},
    }
