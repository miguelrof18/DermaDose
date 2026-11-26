"""Estimated placement of a Slicer surface using anatomical or geometric anchors."""
import numpy as np
from pydicom import dcmread
from trimesh import load_mesh

from .geometry import beam_from_event, table_shift, patient_pose, geometry_preflight


def _coded_values(root, concept):
    values = set()
    for item in getattr(root, "ContentSequence", ()):
        name = getattr(item, "ConceptNameCodeSequence", ())
        value = getattr(item, "ConceptCodeSequence", ())
        if name and value and str(name[0].CodeMeaning).strip().casefold() == concept:
            values.add(str(value[0].CodeMeaning).strip().casefold())
        values.update(_coded_values(item, concept))
    return values


def estimate_placement(path, digest, metadata, dicom, events, equipment):
    if metadata.get("mesh_file") != path.name or metadata.get("mesh_sha256") != digest:
        raise ValueError("El JSON de Slicer no corresponde al STL seleccionado.")
    if metadata.get("coordinate_system") != "RAS" or metadata.get("length_unit") != "mm":
        raise ValueError("La colocación automática requiere un STL documentado en RAS y mm.")
    sr = dcmread(dicom)
    direction = _coded_values(sr, "patient table relationship")
    posture = _coded_values(sr, "patient orientation modifier")
    if direction not in ({"headfirst"}, {"feetfirst"}) or posture not in ({"supine"}, {"prone"}):
        raise ValueError("El RDSR no identifica una postura supina o prona uniforme y el sentido de entrada.")
    regions = _coded_values(sr, "target region")
    anchors = metadata.get("anatomical_anchors_ras_mm", {})
    region = next(iter(regions)) if len(regions) == 1 else None
    mesh = load_mesh(path, process=False)
    bounds = mesh.bounds
    if bounds is None or bounds.shape != (2, 3) or not np.isfinite(bounds).all():
        raise ValueError("La superficie no tiene límites válidos.")
    has_anchor = region is not None and isinstance(anchors, dict) and region in anchors
    anchor = np.asarray(anchors[region], dtype=float) if has_anchor else bounds.mean(axis=0)
    if anchor.shape != (3,) or not np.isfinite(anchor).all():
        raise ValueError("Punto de anclaje inválido.")
    if np.any(anchor < bounds[0]) or np.any(anchor > bounds[1]):
        raise ValueError("El punto de la región diana queda fuera de los límites del STL.")
    # RAS to left, cranial, anterior for head-first supine.
    matrix = np.array([[-1, 0, 0], [0, 0, 1], [0, 1, 0]])
    if direction == {"feetfirst"}:
        matrix = np.diag([-1, -1, 1]) @ matrix
    if posture == {"prone"}:
        matrix = np.diag([-1, 1, -1]) @ matrix
    if not events:
        raise ValueError("No hay eventos para establecer la posición de referencia.")
    for reference_index, event in enumerate(events, 1):
        beam = beam_from_event(event, equipment)
        if abs(beam.direction[2]) >= 0.1:
            break
    else:
        raise ValueError("No hay un haz de referencia suficientemente transversal al plano de mesa.")
    shift = table_shift(event, equipment)
    geometry = equipment.geometry_for(event.plane)
    support_z = float(geometry["table_plane_z_mm"])
    local_anchor = matrix @ anchor
    min_z = float(np.min(mesh.vertices @ matrix[2]))
    offset_z = support_z - min_z
    room_anchor_z = local_anchor[2] + offset_z + shift[2]
    iso = beam.source_mm + float(event.sod_mm) * beam.direction
    depth = (room_anchor_z - beam.source_mm[2]) / beam.direction[2]
    if not 0 < depth < float(event.sid_mm):
        raise ValueError("El anclaje estimado queda fuera del tramo entre foco y detector.")
    target = beam.source_mm + depth * beam.direction
    offset = target - shift - local_anchor
    return dict(mesh_unit="mm", mesh_axes_matrix=matrix.tolist(),
                mesh_origin_mm=offset.tolist(), origin_on_table_mm=[0, 0, 0],
                tilt_x_deg=0, tilt_y_deg=0), dict(
        target_region=region, anchor_ras_mm=anchor.tolist(),
        anchor_source='ct_anatomical' if has_anchor else 'mesh_bounds_center',
        reference_event_index=reference_index, reference_table_shift_mm=shift.tolist(),
        reference_isocenter_mm=iso.tolist(), posture=next(iter(posture)),
        table_relationship=next(iter(direction)),
        assumptions=["superficie apoyada en el plano de mesa configurado",
                     ("anclaje anatómico en el eje central del evento de referencia" if has_anchor else
                      "centro geométrico del volumen delimitador en el eje central del evento de referencia"),
                     "paciente fijo respecto a la mesa"])


def estimate_room_reference(path, metadata, dicom, events, equipment, *, interactive=False, input_fn=input):
    """Estimate only T -> room translation, keeping the manual C -> T registration."""
    if equipment.registered_frames is not None or "ct_table_registration" not in equipment.patient:
        raise ValueError("La estimación requiere un registro TC–tablero provisional.")
    if "table_room_offset_mm" in equipment.geometry:
        return {"method": "configured_room_translation", "offset_room_mm": equipment.geometry["table_room_offset_mm"]}
    sr = dcmread(dicom, stop_before_pixels=True)
    regions = _coded_values(sr, "target region")
    anchors = metadata.get("anatomical_anchors_ras_mm", {})
    if not isinstance(anchors, dict):
        anchors = {}
    region = next(iter(regions)) if len(regions) == 1 else None
    # Use the existing anatomical reference if available; otherwise a clearly
    # labelled geometric proxy. Never invent absolute RAS coordinates or ask
    # for a second manual marker in the ordinary provisional workflow.
    mesh = load_mesh(path, process=False)
    bounds = mesh.bounds
    if bounds is None or bounds.shape != (2, 3) or not np.isfinite(bounds).all():
        raise ValueError("Límites STL inválidos.")
    axes = metadata.get("coordinate_system")
    if axes not in ("RAS", "LPS"):
        raise ValueError("La malla requiere un sistema RAS o LPS documentado.")
    ras_to_native = np.array([1., 1., 1.]) if axes == "RAS" else np.array([-1., -1., 1.])
    thoracic = {"thorax", "thoracic region", "heart", "cardiac", "lung", "lungs", "pulmonary", "chest"}
    key = region if region in anchors else (
        "chest" if "chest" in anchors and (region is None or region in thoracic) else None)
    native = bounds.mean(axis=0)
    ras = native * ras_to_native
    source = "mesh_bounds_center_geometric_only"
    if key is not None:
        try:
            candidate = np.asarray(anchors[key], dtype=float)
        except (TypeError, ValueError):
            candidate = np.empty(0)
        if candidate.shape == (3,) and np.isfinite(candidate).all():
            candidate_native = candidate * ras_to_native
            if np.all(candidate_native >= bounds[0]) and np.all(candidate_native <= bounds[1]):
                ras, native = candidate, candidate_native
                source = "ct_anatomical_" + key
    patient = equipment.patient
    local = np.asarray(patient["mesh_axes_matrix"]) @ native + np.asarray(patient["mesh_origin_mm"])
    # Keep the table height and orientation fixed. Only a horizontal translation is inferred.
    # Prefer a transverse reference, but permit oblique views if the central
    # ray crosses the estimated target height between source and detector.
    for index, event in enumerate(events, 1):
        if event.ka_r_gy is None or event.ka_r_gy <= 0:
            continue
        beam = beam_from_event(event, equipment)
        if abs(beam.direction[2]) < 0.1:
            continue
        shift = table_shift(event, equipment)
        rotation, translation = patient_pose(equipment, shift)
        point = rotation @ local + translation
        depth = (point[2] - beam.source_mm[2]) / beam.direction[2]
        if not 0 < depth < event.sid_mm:
            continue
        central = beam.source_mm + depth * beam.direction
        offset = central - point
        offset[2] = 0.
        equipment.geometry["table_room_offset_mm"] = offset.tolist()
        # The same offset moves the board, its marker and the registered patient.
        evidence = dict(method="target_region_first_suitable_beam", target_source=source,
                        target_region=region, target_ras_mm=ras.tolist(), reference_event_index=index,
                        reference_event_uid=event.uid, offset_room_mm=offset.tolist(),
                        target_before_room_mm=point.tolist(), target_after_room_mm=(point + offset).tolist(),
                        limitation=("Estimated from a coarse target proxy and beam; not measured room calibration. "
                                    "The geometric centre is not a verified organ or irradiation target."))
        patient["ct_table_registration"]["room_reference_estimate"] = evidence
        return evidence
    raise ValueError("No hay un haz de referencia adecuado para estimar la posición del tablero; mida la referencia de sala.")
