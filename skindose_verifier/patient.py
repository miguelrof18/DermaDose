"""Mallas de muestreo locales de la piel; z apunta hacia anterior/techo."""
from dataclasses import dataclass
import math
from pathlib import Path
from functools import cached_property

import numpy as np
import trimesh
from scipy.sparse import csr_matrix

from .domain import InputError
from .mesh_processing import refine_mesh


@dataclass(frozen=True)
class SkinSurface:
    points_mm: np.ndarray  # (N, 3), in patient coordinates
    normals: np.ndarray    # outward unit normals (N, 3)
    area_mm2: np.ndarray   # representative element areas (N,)
    shape: tuple[int, ...]
    axes: tuple[np.ndarray, np.ndarray]  # rows, columns
    axis_labels: tuple[str, str]
    kind: str
    mesh: trimesh.Trimesh | None = None

    @cached_property
    def face_graph(self) -> csr_matrix:
        """Cache adjacency distances once per mesh for neighbourhoods/clusters."""
        if self.mesh is None:
            raise InputError("El grafo de caras requiere malla triangular.")
        a, b = self.mesh.face_adjacency.T
        distance = np.linalg.norm(self.points_mm[a]-self.points_mm[b], axis=1)
        return csr_matrix((np.r_[distance, distance], (np.r_[a, b], np.r_[b, a])),
                          shape=(len(self.points_mm), len(self.points_mm)))


def sample_patient(data: dict) -> SkinSurface:
    kind = data["type"]
    step = float(data["resolution_mm"])
    if not np.isfinite(step) or step <= 0:
        raise InputError("Resolución espacial inválida.")
    if kind == "mesh":
        path = Path(data["mesh_path"])
        if path.suffix.lower() not in (".stl", ".obj", ".ply"):
            raise InputError("Se aceptan mallas triangulares STL, OBJ o PLY.")
        try:
            mesh = trimesh.load_mesh(path, process=True)
        except (OSError, ValueError, TypeError) as exc:
            raise InputError(f"No se pudo leer la malla de piel: {exc}") from exc
        if not isinstance(mesh, trimesh.Trimesh) or len(mesh.faces) < 4:
            raise InputError("Se requiere una única malla de superficie triangular.")
        if not mesh.is_watertight or not mesh.is_winding_consistent:
            raise InputError("Malla abierta o caras con orientación inconsistente; reparar antes de calcular.")
        if len(mesh.split(only_watertight=False)) != 1:
            raise InputError("La malla de piel debe formar una única superficie conectada.")
        unit = data.get("mesh_unit")
        if unit not in ("mm", "cm", "m"):
            raise InputError("Declarar mesh_unit='mm', 'cm' o 'm'; STL no guarda unidades.")
        scale = {"mm": 1., "cm": 10., "m": 1000.}[unit]
        matrix = np.asarray(data.get("mesh_axes_matrix"), dtype=float)
        origin = np.asarray(data.get("mesh_origin_mm"), dtype=float)
        if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)) or not np.allclose(
            matrix.T @ matrix, np.eye(3), atol=1e-7) or not np.isclose(abs(np.linalg.det(matrix)), 1):
            raise InputError("mesh_axes_matrix debe describir una rotación/reflexión ortonormal explícita.")
        if origin.shape != (3,) or not np.all(np.isfinite(origin)):
            raise InputError("mesh_origin_mm inválido.")
        mesh.vertices = (np.asarray(mesh.vertices, dtype=float) * scale) @ matrix.T + origin
        if mesh.volume < 0:
            mesh.invert()
        if not mesh.is_volume or not np.all(np.isfinite(mesh.vertices)):
            raise InputError("Malla sin volumen positivo/normales exteriores válidas.")
        if np.any(mesh.area_faces < 1e-8) or not np.all(np.isfinite(mesh.area_faces)):
            raise InputError("La malla contiene caras degeneradas.")
        if data.get("mesh_max_edge_mm") is not None:
            mesh = refine_mesh(mesh, float(data["mesh_max_edge_mm"]),
                               roi_mm=data.get("mesh_refinement_roi_mm"),
                               max_faces=data.get("mesh_max_faces", 1_000_000))
        points = np.asarray(mesh.triangles_center)
        normals = np.asarray(mesh.face_normals)
        areas = np.asarray(mesh.area_faces)
        if np.any(areas < 1e-8) or not np.all(np.isfinite(areas)):
            raise InputError("La malla contiene caras degeneradas.")
        # One centroid quadrature point per face, including real refined faces.
        return SkinSurface(points, normals, areas, (len(points),),
                           (np.array([]), np.array([])), ("", ""), kind, mesh)
    length = float(data["length_mm"])
    if not np.isfinite(length) or length <= 0:
        raise InputError("Longitud del modelo inválida.")
    ny = max(1, round(length / step))
    dy = length / ny
    yy = -length/2 + (np.arange(ny) + .5)*dy
    if kind == "plane":
        width = float(data["width_mm"])
        if not np.isfinite(width) or width <= 0 or not np.isfinite(float(data["z_mm"])):
            raise InputError("Ancho del plano inválido.")
        nx = max(1, round(width/step))
        dx = width/nx
        xx = -width/2 + (np.arange(nx) + .5)*dx
        x, y = np.meshgrid(xx, yy)
        points = np.column_stack((x.ravel(), y.ravel(), np.full(x.size, float(data["z_mm"]))))
        normal = np.asarray(data["normal"], dtype=float)
        norm = np.linalg.norm(normal)
        if norm == 0 or not np.isfinite(norm):
            raise InputError("Normal del plano inválida.")
        normals = np.broadcast_to(normal / norm, points.shape).copy()
        # Horizontal plane: area is dx*dy; tilted planes require a 2D tangent basis.
        if not np.allclose(np.abs(normal/norm), [0., 0., 1.]):
            raise InputError("El plano discreto requiere normal vertical +/-z.")
        return SkinSurface(points, normals, np.full(points.shape[0], dx*dy),
                           (ny, nx), (yy, xx), ("y longitudinal [mm]", "x lateral [mm]"), kind)
    if kind == "ellipse":
        a, b = float(data["semi_axis_x_mm"]), float(data["semi_axis_z_mm"])
        if not np.isfinite(a+b) or a <= 0 or b <= 0:
            raise InputError("Semiejes de la elipse inválidos.")
        perimeter = math.pi*(3*(a+b) - math.sqrt((3*a+b)*(a+3*b)))
        nt = max(24, round(perimeter/step))
        dt = 2*math.pi/nt
        theta = -math.pi+(np.arange(nt)+.5)*dt
        y, th = np.meshgrid(yy, theta, indexing="ij")
        points = np.column_stack((a*np.cos(th.ravel()), y.ravel(), b*np.sin(th.ravel())))
        normals = np.column_stack((np.cos(th.ravel())/a, np.zeros(th.size), np.sin(th.ravel())/b))
        normals /= np.linalg.norm(normals, axis=1)[:, None]
        d_arc = np.sqrt((a*np.sin(th))**2 + (b*np.cos(th))**2)*dt
        return SkinSurface(points, normals, (d_arc*dy).ravel(), (ny, nt),
                           (yy, theta), ("y longitudinal [mm]", "theta circunferencial [rad]"), kind)
    raise InputError(f"Modelo de piel desconocido: {kind}.")
