"""Estimated placement of a Slicer surface using anatomical or geometric anchors."""
import numpy as np
from pydicom import dcmread
from trimesh import load_mesh

from .geometry import beam_from_event, table_shift


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
