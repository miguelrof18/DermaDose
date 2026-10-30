"""Versioned equipment and factor definitions. No clinical constants in code."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import math
from typing import Any
from functools import cached_property

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .domain import InputError, IrradiationEvent


ALLOWED_AXES = {"kvp", "filtration_mm_cu", "field_area_cm2", "table_angle_deg", "patient_thickness_mm",
                "skin_angle_deg", "field_fraction_x", "field_fraction_y"}


@dataclass(frozen=True)
class Factor:
    name: str
    provenance: str
    identifier: str
    value: float | None
    axes: dict[str, list[float]] | None
    grid: Any
    relative_standard_uncertainty: float | None

    @classmethod
    def from_dict(cls, name: str, data: dict[str, Any]) -> "Factor":
        if not isinstance(data, dict) or not data.get("provenance") or not data.get("identifier"):
            raise InputError(f"Factor {name}: requiere procedencia e identificador.")
        value, axes, grid = data.get("value"), data.get("axes"), data.get("grid")
        if (value is None) == (axes is None):
            raise InputError(f"Factor {name}: definir exactamente value o axes/grid.")
        if value is not None and (not math.isfinite(float(value)) or float(value) <= 0):
            raise InputError(f"Factor {name}: valor no positivo o no finito.")
        if axes is not None:
            if not axes or set(axes) - ALLOWED_AXES or grid is None:
                raise InputError(f"Factor {name}: ejes o tabla inválidos.")
            shape = []
            for key, vals in axes.items():
                if len(vals) < 2 or any(not math.isfinite(float(x)) for x in vals) or any(
                    float(a) >= float(b) for a, b in zip(vals, vals[1:])
                ):
                    raise InputError(f"Factor {name}: eje {key} no creciente.")
                shape.append(len(vals))
            arr = np.asarray(grid, dtype=float)
            if tuple(shape) != arr.shape or not np.all(np.isfinite(arr)) or np.any(arr <= 0):
                raise InputError(f"Factor {name}: dimensiones/valores de tabla inválidos.")
        u = data.get("relative_standard_uncertainty")
        if u is not None and (not math.isfinite(float(u)) or float(u) < 0):
            raise InputError(f"Factor {name}: incertidumbre inválida.")
        return cls(name, str(data["provenance"]), str(data["identifier"]),
                   None if value is None else float(value), axes, grid,
                   None if u is None else float(u))

    @cached_property
    def interpolator(self):
        if self.axes is None:
            return None
        return RegularGridInterpolator(tuple(np.asarray(v, dtype=float) for v in self.axes.values()),
                                       np.asarray(self.grid, dtype=float), bounds_error=True)

    def evaluate(self, **context: float | np.ndarray) -> np.ndarray:
        if self.value is not None:
            return np.asarray(self.value)
        assert self.axes is not None
        try:
            coordinates = np.broadcast_arrays(*(np.asarray(context[k], dtype=float) for k in self.axes))
        except KeyError as exc:
            raise InputError(f"Factor {self.name}: falta parámetro {exc}.") from exc
        if any(not np.all(np.isfinite(c)) for c in coordinates):
            raise InputError(f"Factor {self.name}: argumento no finito.")
        points = np.stack([c.ravel() for c in coordinates], axis=-1)
        try:
            return self.interpolator(points).reshape(coordinates[0].shape)
        except ValueError as exc:
            raise InputError(f"Factor {self.name}: fuera del rango calibrado; no se extrapola.") from exc


@dataclass(frozen=True)
class Equipment:
    raw: dict[str, Any]
    factors: dict[str, Factor]
    version: str
    factors_by_mode: dict[str, dict[str, Factor]]
    surface_factor_provider: Any = None

    @classmethod
    def load(cls, path: str | Path) -> "Equipment":
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.from_dict(data, base_directory=Path(path).resolve().parent)

    @classmethod
    def from_dict(cls, data: dict, *, base_directory: Path = Path(".")) -> "Equipment":
        """Validate a complete configuration, also used for uncertainty draws."""
        from copy import deepcopy
        data = deepcopy(data)
        phantom = data.get("phantom")
        if phantom is not None:
            if not isinstance(phantom, dict) or phantom.get("name") not in ("XCAT50", "ICRP_AM", "ICRP_AF"):
                raise InputError("phantom.name debe ser XCAT50, ICRP_AM o ICRP_AF.")
            if phantom.get("name") != "XCAT50":
                if not phantom.get("mesh_path"):
                    raise InputError(f"El fantoma {phantom['name']} requiere phantom.mesh_path.")
                data["patient"] = {**data.get("patient", {}), "type": "mesh",
                                    "mesh_path": phantom["mesh_path"],
                                    **{k: v for k, v in phantom.items() if k != "name"}}
            data.setdefault("patient", {})["phantom_name"] = phantom.get("name", "XCAT50")
        for key in ("version", "manufacturer", "model", "geometry", "patient", "factors"):
            if key not in data:
                raise InputError(f"Configuración: falta {key}.")
        for flag in ("allow_synthetic_device_override", "allow_historical_factor_fallback", "allow_uniform_sweep"):
            if flag in data and not isinstance(data[flag], bool):
                raise InputError(f"{flag}: requiere booleano JSON.")
        u = data.get("rdsr_relative_standard_uncertainty")
        if u is not None and (not math.isfinite(float(u)) or float(u) < 0):
            raise InputError("Incertidumbre relativa RDSR no finita o negativa.")
        g = data["geometry"]
        required = ("source_to_isocenter_mm", "source_to_detector_mm", "source_to_irp_mm",
                    "field_reference", "primary_sign", "secondary_sign", "rotation_order",
                    "table_coordinates", "table_origin_mm", "table_matrix", "table_plane_z_mm",
                    "table_half_width_mm", "table_half_length_mm")
        if any(k not in g for k in required):
            raise InputError("Geometría incompleta; se requiere caracterización explícita.")
        if g["field_reference"] not in ("detector", "isocenter") or g["rotation_order"] not in ("primary_then_secondary", "secondary_then_primary"):
            raise InputError("Referencia del campo u orden angular inválido.")
        if g["table_coordinates"] not in ("isocenter", "legacy"):
            raise InputError("Sistema de coordenadas de mesa no declarado.")
        M = np.asarray(g["table_matrix"], dtype=float)
        if M.shape != (3, 3) or not np.allclose(M.T @ M, np.eye(3)) or abs(round(np.linalg.det(M))) != 1:
            raise InputError("La matriz de mesa debe ser una permutación/rotación ortonormal.")
        if len(g["table_origin_mm"]) != 3 or not np.all(np.isfinite(g["table_origin_mm"])):
            raise InputError("Origen de mesa inválido.")
        for key in ("source_to_isocenter_mm", "source_to_detector_mm", "source_to_irp_mm"):
            if not math.isfinite(float(g[key])) or float(g[key]) <= 0:
                raise InputError(f"Geometría: {key} inválida.")
        if g["source_to_detector_mm"] <= g["source_to_isocenter_mm"]:
            raise InputError("SID debe superar la distancia fuente-isocentro.")
        if not all(math.isfinite(float(g[k])) for k in ("table_half_width_mm", "table_half_length_mm", "table_plane_z_mm")) or float(g["table_half_width_mm"]) <= 0 or float(g["table_half_length_mm"]) <= 0:
            raise InputError("Dimensiones de mesa inválidas.")
        if data["patient"]["type"] not in ("plane", "ellipse", "mesh"):
            raise InputError("Modelo paciente admitido: plane/ellipse/mesh.")
        if data["patient"]["type"] == "mesh":
            mesh_path = Path(data["patient"].get("mesh_path", ""))
            if not mesh_path.is_absolute():
                mesh_path = base_directory/mesh_path
            if not mesh_path.is_file():
                raise InputError("La malla de piel no existe junto a la configuración.")
            data["patient"]["mesh_path"] = str(mesh_path.resolve())
        factors = {k: Factor.from_dict(k, v) for k, v in data["factors"].items()}
        if "table_total" in factors and ("table_normal" in factors or "table_angular" in factors):
            raise InputError("TAF total excluye TAF normal/angular: doble conteo.")
        modes = {}
        for mode, overrides in data.get("factors_by_mode", {}).items():
            selected = {**factors, **{k: Factor.from_dict(k, v) for k, v in overrides.items()}}
            if "table_total" in selected and ("table_normal" in selected or "table_angular" in selected):
                raise InputError(f"Modo {mode}: TAF total y normal/angular causarían doble conteo.")
            modes[str(mode)] = selected
        provider = None
        if "surface_factors" in data:
            from .surface_factors import make_event_provider
            if data.get("dose_medium") != "water":
                raise InputError("El proveedor superficial agua/aire requiere dose_medium='water'.")
            for block in [factors, *modes.values()]:
                if {"bsf", "f_skin"} & block.keys():
                    raise InputError("surface_factors excluye bsf/f_skin heredados: doble conteo o sustitución ambigua.")
            provider = make_event_provider(data["surface_factors"], base_directory=base_directory)
        return cls(data, factors, str(data["version"]), modes, provider)

    def factors_for(self, ev: IrradiationEvent) -> dict[str, Factor]:
        if not self.factors_by_mode:
            return self.factors
        if not ev.mode or ev.mode not in self.factors_by_mode:
            raise InputError(f"Evento {ev.uid}: modo {ev.mode!r} sin calibración caracterizada.")
        return self.factors_by_mode[ev.mode]

    @property
    def geometry(self) -> dict[str, Any]:
        return self.raw["geometry"]

    def geometry_for(self, plane: str) -> dict[str, Any]:
        planes = self.geometry.get("planes")
        if planes is not None and plane not in planes:
            raise InputError(f"Plano {plane!r} sin geometría de equipo caracterizada.")
        g = {**self.geometry, **(planes[plane] if planes else {})}
        if g["table_coordinates"] not in ("isocenter", "legacy"):
            raise InputError(f"Plano {plane}: sistema de mesa inválido.")
        if g["rotation_order"] not in ("primary_then_secondary", "secondary_then_primary") or g["field_reference"] not in ("detector", "isocenter"):
            raise InputError(f"Plano {plane}: referencia de campo/orden angular inválido.")
        for key in ("isocenter_mm", "table_origin_mm"):
            value = np.asarray(g.get(key, [0, 0, 0]), dtype=float)
            if value.shape != (3,) or not np.all(np.isfinite(value)):
                raise InputError(f"Plano {plane}: {key} inválido.")
        M = np.asarray(g["table_matrix"], dtype=float)
        if M.shape != (3, 3) or not np.all(np.isfinite(M)) or not np.allclose(M.T @ M, np.eye(3), atol=1e-8):
            raise InputError(f"Plano {plane}: matriz de mesa inválida.")
        if any(not math.isfinite(float(g[k])) for k in ("table_plane_z_mm", "table_half_width_mm", "table_half_length_mm")) or min(float(g["table_half_width_mm"]), float(g["table_half_length_mm"])) <= 0:
            raise InputError(f"Plano {plane}: mesa inválida.")
        return g

    @property
    def patient(self) -> dict[str, Any]:
        return self.raw["patient"]

    def require_matching_device(self, manufacturer: str, model: str) -> None:
        if (manufacturer.casefold().strip(), model.casefold().strip()) != (
            str(self.raw["manufacturer"]).casefold().strip(), str(self.raw["model"]).casefold().strip()
        ):
            raise InputError(
                "RDSR leído, pero el fabricante/modelo no coincide con la configuración: "
                f"RDSR={manufacturer!r}/{model!r}; "
                f"configuración={self.raw['manufacturer']!r}/{self.raw['model']!r}. "
                "Se necesita una configuración de geometría y calibración para el equipo del RDSR."
            )

    def validate_event(self, ev: IrradiationEvent, stage: str) -> None:
        fields = ("ka_r_gy", "primary_deg", "secondary_deg", "table_mm", "sid_mm",
                  "sod_mm", "srp_mm", "field_width_mm", "field_height_mm")
        missing = [key for key in fields if getattr(ev, key) is None]
        if missing:
            raise InputError(f"Evento {ev.uid}: faltan {', '.join(missing)}.")
        if ev.table_system != self.geometry_for(ev.plane)["table_coordinates"]:
            raise InputError(f"Evento {ev.uid}: sistema de mesa {ev.table_system} diferente del configurado.")
        if stage not in ("free_air", "calibrated_air", "skin_dose"):
            raise InputError(f"Etapa desconocida: {stage}.")
        factors = self.factors_for(ev) if stage != "free_air" else self.factors
        if stage != "free_air" and ("cf" not in factors or not ({"table_total"} <= factors.keys() or {"table_normal", "table_angular"} <= factors.keys())):
            raise InputError("Etapa calibrada: faltan CF y transmisión de mesa documentados.")
        if stage == "skin_dose" and self.surface_factor_provider is None and not {"bsf", "f_skin"} <= factors.keys():
            raise InputError("Dosis absorbida: faltan BSF y f_skin documentados.")
        if stage == "skin_dose" and self.raw.get("dose_medium") not in ("skin", "water"):
            raise InputError("Dosis absorbida requiere dose_medium='skin' o 'water'.")
        if stage != "free_air" and (ev.kvp is None or ev.filtration_mm_cu is None):
            raise InputError("La etapa calibrada requiere kVp y material/espesor de filtración identificados.")
        if any(not math.isfinite(float(getattr(ev, key))) or float(getattr(ev, key)) <= 0
               for key in ("sid_mm", "sod_mm", "srp_mm", "field_width_mm", "field_height_mm")):
            raise InputError("Distancias y colimación deben ser positivas.")
        if any(not math.isfinite(float(getattr(ev, key))) for key in ("ka_r_gy", "primary_deg", "secondary_deg")):
            raise InputError("Kerma/ángulos del evento no finitos.")
        if float(ev.ka_r_gy) < 0 or (ev.kap_gy_m2 is not None and (not math.isfinite(ev.kap_gy_m2) or ev.kap_gy_m2 < 0)):
            raise InputError("Kerma por evento negativo.")
        if ev.exposure_mas is not None and (not math.isfinite(ev.exposure_mas) or ev.exposure_mas < 0):
            raise InputError("Exposición mAs inválida.")
        if any(x is not None and not math.isfinite(x) for x in (ev.primary_end_deg, ev.secondary_end_deg)):
            raise InputError("Ángulos finales no finitos.")
        if stage != "free_air" and (not math.isfinite(ev.kvp+ev.filtration_mm_cu) or ev.kvp <= 0 or ev.filtration_mm_cu < 0):
            raise InputError("kVp/filtración inválidos.")
