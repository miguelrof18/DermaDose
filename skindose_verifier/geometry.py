"""Geometría C-arm/mesa/paciente. x izquierda paciente; y craneal; z anterior.

Convención matemática de demostración: ángulo primario gira alrededor de +y,
secundario alrededor de +x, haz inicial +z, foco bajo la mesa. Las convenciones
clínicas y el orden se verifican con fantoma y se fijan en cada configuración.
"""
from dataclasses import dataclass
import math

import numpy as np

from .domain import InputError, IrradiationEvent
from .config import Equipment
from .patient import SkinSurface
from .raytrace import first_hits


def rot_x(deg: float) -> np.ndarray:
    a = math.radians(deg)
    c, s = math.cos(a), math.sin(a)
    return np.array([[1., 0., 0.], [0., c, -s], [0., s, c]])


def rot_y(deg: float) -> np.ndarray:
    a = math.radians(deg)
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, 0., s], [0., 1., 0.], [-s, 0., c]])


@dataclass(frozen=True)
class Beam:
    source_mm: np.ndarray
    detector_mm: np.ndarray
    direction: np.ndarray
    transverse_x: np.ndarray
    transverse_y: np.ndarray
    reference_distance_mm: float
    irp_distance_mm: float


@dataclass(frozen=True)
class Projection:
    visible_mask: np.ndarray
    distance_mm: np.ndarray
    table_intersection_mask: np.ndarray
    table_angle_deg: np.ndarray
    field_area_cm2: np.ndarray
    room_points_mm: np.ndarray
    beam: Beam
    patient_translation_mm: np.ndarray
    skin_angle_deg: np.ndarray
    field_fraction_xy: np.ndarray
    statistics: dict


def beam_from_event(ev: IrradiationEvent, equipment: Equipment) -> Beam:
    g = equipment.geometry_for(ev.plane)
    if g["primary_sign"] not in (-1, 1) or g["secondary_sign"] not in (-1, 1):
        raise InputError("Signos angulares deben ser +/-1.")
    a = float(ev.primary_deg) * int(g["primary_sign"])
    b = float(ev.secondary_deg) * int(g["secondary_sign"])
    R = rot_y(a) @ rot_x(b) if g["rotation_order"] == "primary_then_secondary" else rot_x(b) @ rot_y(a)
    direction, x_axis, y_axis = R[:, 2], R[:, 0], R[:, 1]
    sod, sid, srp = float(ev.sod_mm), float(ev.sid_mm), float(ev.srp_mm)
    if sid <= sod or min(sod, sid, srp) <= 0:
        raise InputError("Geometría fuente/isocentro/detector/IRP físicamente inconsistente.")
    iso = np.asarray(g.get("isocenter_mm", [0., 0., 0.]), dtype=float)
    source = iso - direction*sod
    ref = sid if g["field_reference"] == "detector" else sod
    return Beam(source, source+direction*sid, direction, x_axis, y_axis, ref, srp)


def table_shift(ev: IrradiationEvent, equipment: Equipment) -> np.ndarray:
    g = equipment.geometry_for(ev.plane)
    raw = np.asarray(ev.table_mm, dtype=float)
    zero = np.asarray(g["table_origin_mm"], dtype=float)
    if raw.shape != (3,) or not np.all(np.isfinite(raw)):
        raise InputError("Posición de mesa ausente/inválida.")
    return np.asarray(g["table_matrix"], dtype=float) @ (raw - zero)


def beam_view(ev: IrradiationEvent, equipment: Equipment) -> dict:
    """Schematic beam geometry, expressed in the SAME frame as stored dose."""
    b = beam_from_event(ev, equipment)
    rotation = rot_y(float(equipment.patient.get("tilt_y_deg", 0.))) @ rot_x(float(equipment.patient.get("tilt_x_deg", 0.)))
    shift = table_shift(ev, equipment) + np.asarray(equipment.patient.get("origin_on_table_mm", [0, 0, 0]))
    to_patient = lambda points: ((np.asarray(points)-shift) @ rotation).tolist()
    iso = b.source_mm + b.direction*float(ev.sod_mm)
    center = b.source_mm + b.direction*b.reference_distance_mm
    corners = [center + sx*float(ev.field_width_mm)/2*b.transverse_x + sy*float(ev.field_height_mm)/2*b.transverse_y
               for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    arc = [iso + float(ev.sod_mm)*(-np.cos(t)*b.direction + np.sin(t)*b.transverse_x)
           for t in np.linspace(-1.9, 1.9, 61)]
    return {"source_patient_mm": to_patient(b.source_mm), "isocenter_patient_mm": to_patient(iso),
            "detector_patient_mm": to_patient(b.detector_mm), "field_corners_patient_mm": to_patient(corners),
            "primary_deg": ev.primary_deg, "secondary_deg": ev.secondary_deg,
            "table_position_mm": list(ev.table_mm), "gantry_arc_patient_mm": to_patient(arc),
            "diagram_is_schematic": True}


def project(ev: IrradiationEvent, equipment: Equipment, skin: SkinSurface) -> Projection:
    b = beam_from_event(ev, equipment)
    g = equipment.geometry_for(ev.plane)
    shift = table_shift(ev, equipment)
    offset = np.asarray(equipment.patient.get("origin_on_table_mm", [0., 0., 0.]), dtype=float)
    if offset.shape != (3,) or not np.all(np.isfinite(offset)):
        raise InputError("Desplazamiento paciente-mesa inválido.")
    rotation = rot_y(float(equipment.patient.get("tilt_y_deg", 0.))) @ rot_x(float(equipment.patient.get("tilt_x_deg", 0.)))
    points = skin.points_mm @ rotation.T + shift + offset
    normals = skin.normals @ rotation.T
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(normals)):
        raise InputError("Transformación del paciente no finita.")
    source_local = (b.source_mm - shift - offset) @ rotation
    if skin.kind == "ellipse" and abs(source_local[1]) > float(equipment.patient["length_mm"])/2:
        raise InputError("El foco cae fuera de los extremos del cilindro abierto; requieren tapas/malla 3D.")
    delta = points - b.source_mm
    depth = delta @ b.direction
    distance = np.linalg.norm(delta, axis=1)
    front = np.einsum("ij,ij->i", normals, -delta) > 1e-9
    stats = {"rays_traced": 0, "candidate_triangle_pairs": 0,
             "triangles_in_mesh": len(skin.mesh.faces) if skin.mesh is not None else 0,
             "backend": "trimesh_rtree_cpu" if skin.mesh is not None else "analytic_surface"}
    inside = (depth > 0) & (depth <= float(ev.sid_mm)) & front
    half_x = float(ev.field_width_mm)*depth/(2*b.reference_distance_mm)
    half_y = float(ev.field_height_mm)*depth/(2*b.reference_distance_mm)
    inside &= np.abs(delta @ b.transverse_x) <= half_x+1e-9
    inside &= np.abs(delta @ b.transverse_y) <= half_y+1e-9
    if skin.kind == "mesh" and np.any(inside):
        assert skin.mesh is not None
        # Source and triangles are expressed in the unchanged patient frame.
        ids = np.flatnonzero(inside)
        if np.all(source_local > skin.mesh.bounds[0]) and np.all(source_local < skin.mesh.bounds[1]):
            if skin.mesh.contains([source_local])[0]:
                raise InputError("El foco se encuentra dentro de la malla del paciente.")
        first, stats = first_hits(skin.mesh, source_local, skin.points_mm[ids])
        inside[ids[first != ids]] = False
    # Field area belongs to a plane NORMAL to the beam, at axial depth.
    # Euclidean source-point distance would spuriously enlarge it off axis.
    area_cm2 = (float(ev.field_width_mm)*float(ev.field_height_mm)*
                 (depth/b.reference_distance_mm)**2/100.)
    plane_z = float(g["table_plane_z_mm"]) + shift[2]
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (plane_z - b.source_mm[2])/delta[:, 2]
    crossing = np.isfinite(t) & (t > 1e-9) & (t < 1-1e-9)
    t_safe = np.where(crossing, t, 0.)
    xhit = b.source_mm[0]+t_safe*delta[:, 0]
    yhit = b.source_mm[1]+t_safe*delta[:, 1]
    crossing &= np.abs(xhit-shift[0]) <= float(g["table_half_width_mm"])
    crossing &= np.abs(yhit-shift[1]) <= float(g["table_half_length_mm"])
    crossing &= inside
    # Incidence on the tabletop varies across a divergent field; use each ray.
    cos_table = np.clip(np.abs(delta[:, 2])/np.maximum(distance, 1e-12), 0., 1.)
    table_angle = np.degrees(np.arccos(cos_table))
    cosine = np.einsum("ij,ij->i", normals, -delta)/np.maximum(distance, 1e-12)
    skin_angle = np.degrees(np.arccos(np.clip(cosine, -1., 1.)))
    xy = np.column_stack(((delta @ b.transverse_x)/np.maximum(half_x, 1e-12),
                          (delta @ b.transverse_y)/np.maximum(half_y, 1e-12)))
    return Projection(inside, distance, crossing, table_angle, area_cm2, points, b, shift+offset,
                      skin_angle, xy, stats)


def event_slices(ev: IrradiationEvent, equipment: Equipment) -> list[IrradiationEvent]:
    end_a = ev.primary_end_deg if ev.primary_end_deg is not None else ev.primary_deg
    end_b = ev.secondary_end_deg if ev.secondary_end_deg is not None else ev.secondary_deg
    if math.isclose(float(end_a), float(ev.primary_deg), abs_tol=1e-9) and math.isclose(
        float(end_b), float(ev.secondary_deg), abs_tol=1e-9):
        return [ev]
    if not equipment.raw.get("allow_uniform_sweep", False):
        raise InputError(f"Evento {ev.uid}: barrido angular sin perfil temporal, bloqueado.")
    from dataclasses import replace
    n = int(equipment.raw.get("sweep_samples", 0))
    if n < 2:
        raise InputError("Barrido: sweep_samples debe ser >=2.")
    # Explicit research-only assumption: dose rate constant over sweep.
    return [replace(ev.scaled(1/n), primary_deg=float(ev.primary_deg)+(float(end_a)-float(ev.primary_deg))*(k+.5)/n,
                    secondary_deg=float(ev.secondary_deg)+(float(end_b)-float(ev.secondary_deg))*(k+.5)/n,
                    primary_end_deg=None, secondary_end_deg=None)
            for k in range(n)]
