"""Rigid CT (C, LPS) → tabletop (T, head edge) → isocentre (I), in mm.

Frames (docs/sistemas_de_referencia.md):
  C  physical DICOM coordinates of the CT, LPS: x left, y posterior, z head.
  T  tabletop: origin at the centre of the head edge on the upper face;
     x left of the board, y down, z towards the head edge.
  I  isocentre of the C-arm; axes parallel to T with a level, unrotated board.

  p_T = R_T<-C p_C + t_T<-C            (patient placement, per case)
  t_k = t_0 + S (m_k - m_0)            (O_T in I for event k)
  p_I(k) = p_T + t_k                   (R_I<-T = Id: level, unrotated board)

Every link carries its evidence: explicit/measured or estimated. Estimates are
never silent: they are stored in the configuration with method and hypothesis.
"""
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from .domain import InputError


# I = (left, down, head); the existing ray engine uses (left, head, up).
ROOM_FROM_I = np.array([[1., 0., 0.], [0., 0., 1.], [0., -1., 0.]])
RAS_TO_LPS = np.diag([-1., -1., 1.])


def _settings(value, directory, name):
    if isinstance(value, str):
        path = Path(directory) / value
        try:
            value = json.loads(path.read_text(encoding="utf-8-sig"))
        except (OSError, ValueError) as exc:
            raise InputError(f"No se pudo leer {name} en {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise InputError(f"{name} debe ser un objeto JSON o la ruta de su archivo.")
    return value


def _array(value, shape, name):
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise InputError(f"{name}: se requieren valores numéricos.") from exc
    if array.shape != shape or not np.isfinite(array).all():
        raise InputError(f"Falta completar {name} con valores finitos de forma {shape}.")
    return array


def rigid(rotation, translation):
    matrix = np.eye(4)
    matrix[:3, :3] = rotation
    matrix[:3, 3] = translation
    return matrix


def configure_table(geometry, base_directory):
    """Use measured dimensions and the detector directions specified for this equipment."""
    table = _settings(geometry["table_calibration"], base_directory, "table_calibration")
    if not table.get("id"):
        raise InputError("La calibración de mesa requiere un identificador.")
    length, width, thickness = _array([table.get(k) for k in ("length_mm", "width_mm", "thickness_mm")],
                                      (3,), "dimensiones del tablero (mm)")
    if min(length, width, thickness) <= 0:
        raise InputError("Las dimensiones del tablero deben ser positivas.")
    table.update(length_mm=float(length), width_mm=float(width), thickness_mm=float(thickness))
    # Signs refer to rot_y and rot_x in the engine, not to DICOM axis names.
    try:
        geometry["primary_sign"] = {"left": 1, "right": -1}[table["positive_primary_detector"]]
        geometry["secondary_sign"] = {"cranial": -1, "caudal": 1}[table["positive_secondary_detector"]]
    except (KeyError, TypeError) as exc:
        raise InputError("Indique los sentidos positivos del detector para ambos ángulos.") from exc
    geometry.update(table_calibration=table, table_half_length_mm=length / 2, table_half_width_mm=width / 2)


READING_NAMES = ("longitudinal (DCM 113751)", "lateral (DCM 113752)", "altura (DCM 113753)")


def _matrix_with_unknown(value, name):
    """3x3 matrix where JSON null marks a component not yet characterized (NaN)."""
    if not isinstance(value, (list, tuple)) or len(value) != 3 or any(
            not isinstance(row, (list, tuple)) or len(row) != 3 for row in value):
        raise InputError(f"Falta completar {name} como matriz 3 × 3.")
    try:
        matrix = np.array([[np.nan if v is None else float(v) for v in row] for row in value])
    except (TypeError, ValueError) as exc:
        raise InputError(f"{name}: se requieren valores numéricos o null.") from exc
    if np.any(np.isinf(matrix)):
        raise InputError(f"{name}: valores no finitos.")
    return matrix


def unknown_motion_columns(matrix):
    """Readings whose physical direction in I is not characterized."""
    return [j for j in range(3) if not np.all(np.isfinite(matrix[:, j]))]


def motion_matrix(table):
    """Columns correspond to DCM 113751, 113752 and 113753 respectively.

    A column containing null is an uncharacterized reading. It is accepted here
    and rejected later only if that reading changes with respect to m_0.
    """
    trials = table.get("axis_trials")
    if trials is None:
        matrix = _matrix_with_unknown(table.get("motion_matrix_i_per_mm"), "S (motion_matrix_i_per_mm)")
        known = [j for j in range(3) if j not in unknown_motion_columns(matrix)]
        if not known:
            raise InputError("S no caracteriza ninguna lectura de mesa.")
        if np.linalg.matrix_rank(matrix[:, known]) != len(known):
            raise InputError("S debe caracterizar desplazamientos independientes.")
        return matrix
    else:
        if not isinstance(trials, list) or len(trials) != 3:
            raise InputError("S requiere tres movimientos independientes, uno por lectura.")
        columns = {}
        for trial in trials:
            if not isinstance(trial, dict):
                raise InputError("Cada movimiento de mesa debe contener las diferencias observadas.")
            change = _array(trial.get("delta_long_lat_height_mm"), (3,), "cambio de lecturas")
            moved = np.flatnonzero(~np.isclose(change, 0., atol=1e-6))
            if len(moved) != 1 or moved[0] in columns:
                raise InputError("Mueva un solo eje distinto en cada medida de S.")
            index = moved[0]
            physical = _array(trial.get("delta_isocenter_mm"), (3,), "movimiento físico en I")
            columns[index] = physical / change[index]
        matrix = np.column_stack([columns[j] for j in range(3)])
    if np.linalg.matrix_rank(matrix) != 3:
        raise InputError("S debe caracterizar tres desplazamientos independientes.")
    return matrix


def ct_to_table(registration):
    """Fit labelled fiducials, or place a single anchor with a supplied rotation."""
    if "fiducials_ct_lps_mm" in registration or "fiducials_table_mm" in registration:
        if any(key in registration for key in ("anchor_ct_lps_mm", "anchor_table_mm", "rotation_table_from_ct")):
            raise InputError("Seleccione un único método de registro TC–tablero.")
        markers = registration.get("fiducials_ct_lps_mm")
        if not isinstance(markers, list) or len(markers) < 3:
            raise InputError("El registro con marcas requiere tres o más correspondencias 3D.")
        ct = _array(markers, (len(markers), 3), "marcas del TC en LPS")
        table = _array(registration.get("fiducials_table_mm"), ct.shape, "marcas del tablero")
        x, y = ct - ct.mean(axis=0), table - table.mean(axis=0)
        covariance = x.T @ y
        if np.linalg.matrix_rank(covariance) < 2:
            raise InputError("Las marcas deben formar al menos un triángulo no degenerado.")
        u, _, vt = np.linalg.svd(covariance)
        correction = np.diag([1., 1., np.sign(np.linalg.det(vt.T @ u.T))])
        rotation = vt.T @ correction @ u.T
        offset = table.mean(axis=0) - rotation @ ct.mean(axis=0)
        rms = float(np.sqrt(np.mean(np.sum((ct @ rotation.T + offset - table)**2, axis=1))))
        return rigid(rotation, offset), {"method": "rigid_fiducials", "count": len(ct), "fit_rms_mm": rms}
    rotation = _array(registration.get("rotation_table_from_ct"), (3, 3), "rotación TC a tablero")
    if not (np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-6, rtol=0)
            and np.isclose(np.linalg.det(rotation), 1., atol=1e-6, rtol=0)):
        raise InputError("La orientación TC–tablero debe ser una rotación propia.")
    a_ct = _array(registration.get("anchor_ct_lps_mm"), (3,), "anclaje TC en LPS y mm")
    a_table = _array(registration.get("anchor_table_mm"), (3,), "anclaje del tablero en mm")
    return rigid(rotation, a_table - rotation @ a_ct), {"method": "anchor_and_rotation"}


def marker_in_isocenter(reference):
    """Localize the calibration marker from calibrated pinhole views in the same pose."""
    views = reference.get("calibrated_views")
    if views is None:
        return _array(reference.get("marker_isocenter_mm"), (3,), "el mismo marcador en I (mm)"), {}
    if reference.get("marker_isocenter_mm") is not None or not isinstance(views, list) or len(views) < 2:
        raise InputError("Use un punto 3D medido o al menos dos proyecciones calibradas del marcador.")
    projectors, projected_sources, cameras = [], [], []
    for view in views:
        if not isinstance(view, dict):
            raise InputError("Cada vista requiere la matriz de proyección y el píxel del marcador.")
        p = _array(view.get("projection_px_from_i_mm"), (3, 4), "proyección calibrada 3×4")
        uv = _array(view.get("pixel_uv"), (2,), "píxel del marcador (u, v)")
        try:
            source = np.linalg.solve(p[:, :3], -p[:, 3])
            ray = np.linalg.solve(p[:, :3], [*uv, 1.])
        except np.linalg.LinAlgError as exc:
            raise InputError("La matriz de proyección debe representar una cámara de foco finito.") from exc
        direction = ray / np.linalg.norm(ray)
        normal_plane = np.eye(3) - np.outer(direction, direction)
        projectors.append(normal_plane)
        projected_sources.append(normal_plane @ source)
        cameras.append((p, uv))
    # Least-squares point closest to the back-projected lines.
    q, _, rank, _ = np.linalg.lstsq(np.vstack(projectors), np.concatenate(projected_sources), rcond=1e-8)
    if rank != 3:
        raise InputError("Las vistas del marcador carecen de separación angular suficiente.")
    squared_errors = []
    for p, uv in cameras:
        h = p @ np.r_[q, 1.]
        if h[2] <= np.finfo(float).eps * np.linalg.norm(p[2]):
            raise InputError("El marcador debe quedar delante de cada foco (profundidad positiva).")
        squared_errors.append(np.sum((h[:2] / h[2] - uv)**2))
    return q, {"marker_views": len(views), "marker_reprojection_rms_px": float(np.sqrt(np.mean(squared_errors)))}


REFERENCE_KEYS = ("readings_long_lat_height_mm", "marker_table_mm", "marker_isocenter_mm",
                  "calibrated_views", "origin_isocenter_mm")


def room_reference_state(table):
    """'missing', 'measured' (marker in T and I) or 'estimated' (O_T in I given directly)."""
    if not isinstance(table, dict):
        return "missing"
    reference = table.get("reference")
    if reference is None or (isinstance(reference, dict) and all(reference.get(k) is None for k in REFERENCE_KEYS)):
        return "missing"
    if not isinstance(reference, dict):
        raise InputError("table_calibration.reference debe ser un objeto JSON.")
    if reference.get("origin_isocenter_mm") is not None:
        if any(reference.get(k) is not None for k in ("marker_table_mm", "marker_isocenter_mm", "calibrated_views")):
            raise InputError("La referencia de sala es medida (marcador) o estimada (origin_isocenter_mm), no ambas.")
        if not str(reference.get("status", "")).startswith("estimated"):
            raise InputError("origin_isocenter_mm requiere status 'estimated_…' y su método en 'estimate'.")
        return "estimated"
    return "measured"


def table_readings(event):
    """(longitudinal, lateral, height) in mm from an event in the legacy reading system."""
    if event.table_system != "legacy":
        raise InputError(f"Evento {event.uid}: faltan las lecturas DCM 113751–113753 de la mesa calibrada.")
    # IrradiationEvent retains the historical (lateral, longitudinal, height) order.
    lateral, longitudinal, height = _array(event.table_mm, (3,), "posición de mesa del evento")
    return np.array([longitudinal, lateral, height])


def validate_slice_geometry(ct):
    """Verify the entire DICOM stack against the voxel-to-LPS affine, independent of UID."""
    slices = ct.get("dicom_image_geometry")
    if not isinstance(slices, list) or len(slices) < 3 or not all(isinstance(v, dict) for v in slices):
        raise InputError("Se requieren al menos tres cortes con geometría DICOM verificable.")
    if ct.get("unavailable_dicom_instances", 0):
        raise InputError("Faltan imágenes DICOM del volumen.")
    n = len(slices)
    a = _array(ct.get("ijk_to_lps_mm"), (4, 4), "matriz IJK a LPS")
    if not np.allclose(a[3], [0, 0, 0, 1]) or abs(np.linalg.det(a[:3, :3])) < 1e-10:
        raise InputError("Matriz IJK a LPS inválida.")
    iop = _array([v.get("image_orientation_patient") for v in slices], (n, 6), "IOP de los cortes")
    ipp = _array([v.get("image_position_patient_mm") for v in slices], (n, 3), "IPP de los cortes")
    row, col = iop[:, :3], iop[:, 3:]
    if (not np.allclose(np.linalg.norm(row, axis=1), 1., atol=1e-4, rtol=0)
            or not np.allclose(np.linalg.norm(col, axis=1), 1., atol=1e-4, rtol=0)
            or np.max(np.abs(np.sum(row * col, axis=1))) > 1e-4
            or np.max(np.linalg.norm(row - row[0], axis=1)) > 1e-3
            or np.max(np.linalg.norm(col - col[0], axis=1)) > 1e-3):
        raise InputError("Orientación de los cortes no uniforme u ortogonal.")
    i, j, k = a[:3, :3].T
    si, sj = np.linalg.norm(i), np.linalg.norm(j)
    if np.linalg.norm(i / si - row[0]) > 1e-3 or np.linalg.norm(j / sj - col[0]) > 1e-3:
        raise InputError("La orientación DICOM no coincide con la matriz del volumen Slicer.")
    normal = np.cross(row[0], col[0])
    step = float(k @ normal)
    if abs(step) < 1e-5:
        raise InputError("Separación de cortes inválida.")
    indices = (ipp - a[:3, 3]) @ normal / step
    nearest = np.rint(indices).astype(int)
    if not np.array_equal(np.sort(nearest), np.arange(n)):
        raise InputError("Cortes duplicados, ausentes o fuera de la rejilla del volumen.")
    error = np.linalg.norm(ipp - (a[:3, 3] + nearest[:, None] * k), axis=1)
    tolerance = min(0.5, max(0.05, 0.1 * abs(step)))
    if float(error.max()) > tolerance:
        raise InputError("Posiciones de los cortes incompatibles con la rejilla afín del TC.")
    spacings = [v.get("pixel_spacing_mm") for v in slices]
    spacing_checked = any(v is not None for v in spacings)
    if spacing_checked:
        spacing = _array(spacings, (n, 2), "PixelSpacing de los cortes")
        if np.any(spacing <= 0) or not np.allclose(spacing, [sj, si], atol=0.01, rtol=0.001):
            raise InputError("PixelSpacing no uniforme o incompatible con Slicer.")
    dims = ct.get("volume_dimensions_ijk")
    if dims is not None:
        dims = _array(dims, (3,), "dimensiones IJK")
        if np.any(dims <= 0) or not np.array_equal(dims, np.rint(dims)) or dims[2] != n:
            raise InputError("Número de cortes incompatible con el volumen.")
    sizes_checked = any(v.get("rows") is not None or v.get("columns") is not None for v in slices)
    if sizes_checked:
        sizes = _array([[v.get("columns"), v.get("rows")] for v in slices], (n, 2), "dimensiones DICOM")
        if np.any(sizes <= 0) or not np.array_equal(sizes, np.rint(sizes)) or not np.all(sizes == sizes[0]):
            raise InputError("Dimensiones DICOM no uniformes.")
        if dims is not None and not np.array_equal(sizes[0], dims[:2]):
            raise InputError("Dimensiones DICOM incompatibles con Slicer.")
    return dict(status="slice_grid_consistent", slice_count=n, slice_spacing_mm=abs(step),
                max_position_residual_mm=float(error.max()), position_tolerance_mm=tolerance,
                pixel_spacing_checked=spacing_checked, image_dimensions_checked=sizes_checked,
                volume_dimensions_checked=dims is not None,
                limitation="Coherencia geométrica; no acredita identidad del paciente ni anonimización.")


POSE_ROTATIONS = {"HFS": np.eye(3), "FFS": np.diag([-1., 1., -1.]),
                  "HFP": np.diag([-1., -1., 1.]), "FFP": np.diag([1., -1., -1.])}


def intervention_pose(events):
    values = [e.original.get("intervention_patient_position") for e in events]
    known = {v for v in values if v is not None}
    if len(known) > 1:
        raise InputError("El RDSR contiene posiciones distintas; un único registro rígido no es válido.")
    pose = next(iter(known)) if known else None
    if pose and pose not in POSE_ROTATIONS:
        raise InputError("Postura del RDSR no soportada por este registro.")
    return dict(position=pose, event_count=len(values), known_event_count=sum(v is not None for v in values),
                status="uniform" if pose and all(v == pose for v in values) else "partial" if pose else "missing")


def check_registration_pose(registration, events):
    pose = intervention_pose(events)
    transform, _ = ct_to_table(registration)
    r = transform[:3, :3]
    if pose["position"]:
        expected = POSE_ROTATIONS[pose["position"]]
        if np.dot(r[:, 2], expected[:, 2]) <= 0 or np.dot(r[:, 1], expected[:, 1]) <= 0:
            raise InputError("La orientación del registro TC–tablero contradice la postura del RDSR.")
        pose["registration_comparison"] = "compatible_nominal_axes"
    else:
        pose["registration_comparison"] = "manual_rotation_retained_rdsr_missing"
    return pose


MIN_SHOULDER_TO_HEAD_EDGE_MM = 150.


def check_shoulder_anchor_distance(registration) -> None:
    """A shoulder anchor a few centimetres from the head edge leaves the head off the board."""
    shoulder = (registration.get("anchor_definition") == "posterior_shoulder_midline"
                or registration.get("orientation_basis") == "operator_confirmed_aligned_head_first_supine")
    a_table = registration.get("anchor_table_mm")
    if shoulder and a_table is not None and len(a_table) == 3:
        distance = -float(a_table[2])
        if distance < MIN_SHOULDER_TO_HEAD_EDGE_MM:
            raise InputError(
                f"registro_tc.json sitúa los hombros a {distance:g} mm del borde de cabecera del tablero; con la "
                f"cabeza sobre la mesa se esperan al menos {MIN_SHOULDER_TO_HEAD_EDGE_MM:g} mm (típico, 300–450 mm). "
                "Las medidas del tablero van en mm: 40 cm son 400 mm.")


def resolve_registration_pose(registration, events, *, interactive=False, input_fn=input):
    registration = dict(registration)
    check_shoulder_anchor_distance(registration)
    pose = intervention_pose(events)
    try:
        check_registration_pose(registration, events)
    except InputError:
        expected = POSE_ROTATIONS.get(pose["position"])
        nominal = registration.get("orientation_basis") == "operator_confirmed_aligned_head_first_supine"
        if not (interactive and pose["status"] == "uniform" and nominal
                and "anchor_ct_lps_mm" in registration and "anchor_table_mm" in registration
                and expected is not None and np.allclose(registration.get("rotation_table_from_ct"), np.eye(3))):
            raise
        print("Registro manual HFS; RDSR " + pose["position"] + ".")
        answer = input_fn("¿Confirma la postura del intervencionismo del RDSR conservando ambas marcas? [s/N]: ").strip().casefold()
        if answer not in ("s", "si", "sí"):
            raise
        registration["rotation_table_from_ct"] = expected.tolist()
        registration["orientation_basis"] = "operator_confirmed_rdsr_nominal_pose"
    registration["intervention_pose_audit"] = check_registration_pose(registration, events)
    return registration


def validated_ct_to_table(table, patient, base_directory):
    """Validate the patient registration independently of the room calibration."""
    case = _settings(patient["ct_registration"], base_directory, "ct_registration")
    transform, evidence = ct_to_table(case)
    mesh = Path(patient["mesh_path"])
    meta = _settings(str(mesh.with_suffix(".json")), base_directory, "metadatos de la malla TC")
    mesh_hash = sha256()
    with mesh.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            mesh_hash.update(chunk)
    digest = mesh_hash.hexdigest()
    if meta.get("mesh_file") != mesh.name or meta.get("mesh_sha256") != digest:
        raise InputError("Los metadatos de la malla TC no corresponden al STL seleccionado.")
    position_evidence = meta.get("patient_position_evidence", {})
    ct = position_evidence.get("ct", {}) if isinstance(position_evidence, dict) else None
    if not isinstance(ct, dict):
        raise InputError("Los metadatos del marco de referencia del TC son inválidos.")
    uid = ct.get("frame_of_reference_uid") or meta.get("frame_of_reference_uid")
    frame_uids = ct.get("frame_of_reference_uids")
    if frame_uids is not None and (not isinstance(frame_uids, list) or
                                   any(not isinstance(v, str) or not v.strip() for v in frame_uids)):
        raise InputError("Lista de UID del marco de referencia inválida.")
    known_multiple = frame_uids is not None and len(set(frame_uids)) > 1
    if bool(ct.get("mixed_frames")) or known_multiple:
        if case.get("mesh_sha256") != digest:
            raise InputError("Los UID ambiguos requieren vincular el registro al SHA-256 de la malla.")
        if case.get("frame_of_reference_uid"):
            raise InputError("No se puede asociar un UID único a una serie con UID ambiguos.")
        evidence["ct_slice_geometry"] = validate_slice_geometry(ct)
        evidence["frame_of_reference_status"] = "uid_ambiguous_slice_grid_consistent"
        evidence["distinct_frame_uid_count"] = len(set(frame_uids)) if frame_uids is not None else None
    elif ct.get("dicom_image_geometry"):
        evidence["ct_slice_geometry"] = validate_slice_geometry(ct)
        evidence["frame_of_reference_status"] = ("incomplete_dicom_uids" if ct.get("missing_frame_of_reference_uids")
                                                 else "slice_grid_consistent")
    elif ct.get("missing_frame_of_reference_uids"):
        evidence["frame_of_reference_status"] = "incomplete_dicom_uids"
    if not case.get("mesh_sha256") and not case.get("frame_of_reference_uid"):
        raise InputError("Asocie el registro al TC con mesh_sha256 o frame_of_reference_uid.")
    if "mesh_sha256" in case and case["mesh_sha256"] != digest:
        raise InputError("El registro corresponde a otra malla.")
    if "frame_of_reference_uid" in case and case["frame_of_reference_uid"] != uid:
        raise InputError("El Frame of Reference UID del registro y el de la malla no coinciden.")
    axes = meta.get("coordinate_system")
    if axes not in ("LPS", "RAS") or meta.get("length_unit") != "mm":
        raise InputError("La malla TC requiere coordenadas RAS o LPS documentadas en mm.")
    if any(patient.get(k, 0) for k in ("tilt_x_deg", "tilt_y_deg")):
        raise InputError("Declare la orientación completa en el registro TC–tablero.")
    patient.update(ct_registration=case, mesh_unit="mm", mesh_origin_mm=[0., 0., 0.],
                   mesh_axes_matrix=(RAS_TO_LPS if axes == "RAS" else np.eye(3)).tolist())
    status = str(case.get("status", "explicit"))
    if status != "explicit" and not status.startswith("estimated"):
        raise InputError("ct_registration.status debe ser 'explicit' o 'estimated_…'.")
    evidence.update(ct_to_table_status=status, ct_to_table_estimate=case.get("estimate"),
                    M_T_from_C=transform.tolist())
    evidence.update(mesh_input_axes=axes, surface_frame="CT_LPS_mm", mesh_sha256=digest,
                    frame_of_reference_uid=uid, table_calibration_id=table["id"],
                    ct_acquisition_position=ct.get("patient_position"),
                    intervention_pose=case.get("intervention_pose_audit"))
    return transform, evidence


@dataclass(frozen=True)
class RegisteredFrames:
    table_from_ct: np.ndarray        # M_T<-C (4 x 4)
    reference_readings: np.ndarray   # m_0, (long, lat, height) mm
    reference_origin_i: np.ndarray   # t_0, O_T in I at m_0, mm
    displacement_matrix: np.ndarray  # S, mm in I per mm of reading; NaN = unknown
    evidence: dict

    @classmethod
    def from_inputs(cls, geometry, patient, base_directory):
        table = geometry.get("table_calibration")
        if not isinstance(table, dict) or patient["type"] != "mesh":
            raise InputError("El registro TC requiere una malla y geometry.table_calibration.")
        if geometry.get("table_coordinates") != "legacy":
            raise InputError("Esta calibración usa las lecturas DCM 113751–113753 (table_coordinates=legacy).")
        state = room_reference_state(table)
        if state == "missing":
            raise InputError(
                "Falta la referencia tablero–isocentro (t0, m0). Mídala con un marcador "
                "(table_calibration.reference) o deje que el lanzador la estime: "
                "python -m skindose_verifier.example_runner CASO.")
        ref = table["reference"]
        m0 = _array(ref.get("readings_long_lat_height_mm"), (3,), "m0 longitudinal/lateral/altura (mm)")
        s = motion_matrix(table)
        if state == "measured":
            q_t = _array(ref.get("marker_table_mm"), (3,), "marcador en T (mm)")
            q_i, reference_evidence = marker_in_isocenter(ref)
            origin = q_i - q_t
            reference_evidence.update(room_reference_status="measured_marker",
                                      marker_table_mm=q_t.tolist(), marker_isocenter_mm=q_i.tolist())
        else:
            origin = _array(ref.get("origin_isocenter_mm"), (3,), "origen del tablero en I (mm)")
            reference_evidence = dict(room_reference_status=str(ref["status"]),
                                      room_reference_estimate=ref.get("estimate"))
        transform, evidence = validated_ct_to_table(table, patient, base_directory)
        status = table.get("motion_matrix_status")
        evidence.update(origin_reference_i_mm=origin.tolist(), readings_reference_mm=m0.tolist(),
                        motion_matrix_i_per_mm=[[None if not np.isfinite(v) else float(v) for v in row] for row in s],
                        motion_matrix_unknown_readings=[READING_NAMES[j] for j in unknown_motion_columns(s)],
                        motion_matrix_status=status, **reference_evidence)
        return cls(transform, m0, origin, s, evidence)

    def table_origin_i(self, event):
        """t_k = t_0 + S (m_k - m_0): position of O_T in I during event k."""
        delta = table_readings(event) - self.reference_readings
        for j in unknown_motion_columns(self.displacement_matrix):
            if abs(delta[j]) > 1e-6:
                raise InputError(
                    f"Evento {event.uid}: la lectura {READING_NAMES[j]} cambia {delta[j]:+.1f} mm respecto a la "
                    "referencia, pero su sentido en el angiógrafo no está caracterizado (columna null de "
                    "motion_matrix_i_per_mm). Obtenga S con las tomas 4–6 de «Calibración geométrica» "
                    "(python -m skindose_verifier.table_calibration) o declare el sentido comprobado.")
        known = np.where(np.isfinite(self.displacement_matrix), self.displacement_matrix, 0.)
        return self.reference_origin_i + known @ delta
