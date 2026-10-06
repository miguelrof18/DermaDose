"""Objetos físicos en mm, grados, Gy y Gy m²; valores ausentes son None."""
from dataclasses import dataclass, field, replace
from typing import Any


@dataclass(frozen=True)
class IrradiationEvent:
    uid: str
    plane: str
    mode: str | None = None
    ka_r_gy: float | None = None
    kap_gy_m2: float | None = None
    kvp: float | None = None
    filtration_mm_cu: float | None = None
    exposure_mas: float | None = None
    primary_deg: float | None = None
    secondary_deg: float | None = None
    primary_end_deg: float | None = None
    secondary_end_deg: float | None = None
    table_mm: tuple[float, float, float] | None = None
    table_system: str | None = None
    sid_mm: float | None = None
    sod_mm: float | None = None
    srp_mm: float | None = None
    field_width_mm: float | None = None
    field_height_mm: float | None = None
    timestamp: str | None = None
    quality: dict[str, str] = field(default_factory=dict)
    original: dict[str, Any] = field(default_factory=dict)

    def scaled(self, fraction: float) -> "IrradiationEvent":
        """Assign an exact fraction of integrated event kerma to a sweep sample."""
        return replace(self,
            ka_r_gy=None if self.ka_r_gy is None else self.ka_r_gy * fraction,
            kap_gy_m2=None if self.kap_gy_m2 is None else self.kap_gy_m2 * fraction,
            exposure_mas=None if self.exposure_mas is None else self.exposure_mas * fraction)


@dataclass(frozen=True)
class RDSRData:
    events: tuple[IrradiationEvent, ...]
    manufacturer: str
    model: str
    sha256: str
    accumulated_ka_r_gy: float | None = None
    accumulated_kap_gy_m2: float | None = None
    warnings: tuple[str, ...] = ()


class InputError(ValueError):
    """A critical input is missing, inconsistent, or not characterized."""
