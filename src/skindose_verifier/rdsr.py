"""Lectura de contenido SR por (CodingSchemeDesignator, CodeValue), nunca por texto.

Los campos SR están anidados en ContentSequence (0040,A730). Los códigos DCM
siguen TID 10003/10003B/10003C y CID 10008, DICOM PS3.16. No se confunden
los valores por evento con los valores acumulados del estudio.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import replace
from hashlib import sha256
from pathlib import Path
import math
from typing import Iterable

import pydicom
from pydicom.dataset import Dataset

from .domain import InputError, IrradiationEvent, RDSRData

XRAY_DOSE_SR_UID = "1.2.840.10008.5.1.4.1.1.88.67"
EVENT = "113706"
CODE = {
    "plane": "113764", "uid": "113769", "mode": "113721", "timestamp": "111526",
    "ka_r_gy": "113738", "kap_gy_m2": "122130", "kvp": "113733",
    "exposure_mas": "113736", "exposure_time": "113824", "current_ma": "113734",
    "filtration_mm_cu": "113758", "filter_material": "113757",
    "primary_deg": "112011", "secondary_deg": "112012",
    "primary_end_deg": "113739", "secondary_end_deg": "113740",
    "table_x_iso": "128766", "table_y_iso": "128767", "table_z_iso": "128768",
    "table_long": "113751", "table_lat": "113752", "table_height": "113753",
    "sid_mm": "113750", "sod_mm": "113748", "srp_mm": "113737",
    "field_width_mm": "113789", "field_height_mm": "113788",
    "field_area_m2": "113790",
}
UNITS = {
    "ka_r_gy": {"Gy": 1., "mGy": 1e-3, "uGy": 1e-6},
    "kap_gy_m2": {"Gy.m2": 1., "Gy.cm2": 1e-4, "cGy.cm2": 1e-6, "uGy.m2": 1e-6},
    "kvp": {"kV": 1.}, "exposure_mas": {"uA.s": 1e-3, "mA.s": 1.},
    "exposure_time": {"ms": 1., "s": 1000.}, "current_ma": {"mA": 1.},
    "filtration_mm_cu": {"mm": 1.}, "field_area_m2": {"m2": 1., "cm2": 1e-4},
    **{name: {"mm": 1., "cm": 10.} for name in (
        "table_x_iso", "table_y_iso", "table_z_iso", "table_long", "table_lat",
        "table_height", "sid_mm", "sod_mm", "srp_mm", "field_width_mm", "field_height_mm")},
    **{name: {"deg": 1.} for name in (
        "primary_deg", "secondary_deg", "primary_end_deg", "secondary_end_deg")},
}


def concept(item: Dataset) -> tuple[str, str] | None:
    seq = getattr(item, "ConceptNameCodeSequence", None)
    if not seq:
        return None
    return str(seq[0].CodingSchemeDesignator), str(seq[0].CodeValue)


def descendants(item: Dataset):
    for child in getattr(item, "ContentSequence", []):
        yield child
        yield from descendants(child)


def nodes(item: Dataset, code: str) -> list[Dataset]:
    return [node for node in descendants(item) if concept(node) == ("DCM", code)]


def numeric(node: Dataset, field: str) -> tuple[float, str]:
    if getattr(node, "ValueType", None) != "NUM":
        raise InputError(f"{field}: se esperaba un valor NUM.")
    seq = getattr(node, "MeasuredValueSequence", None)
    if not seq or len(seq) != 1 or not getattr(seq[0], "MeasurementUnitsCodeSequence", None):
        raise InputError(f"{field}: NUM sin valor o unidad explícitos.")
    measured = seq[0]
    unit_code = measured.MeasurementUnitsCodeSequence[0]
    if str(unit_code.CodingSchemeDesignator) != "UCUM":
        raise InputError(f"{field}: unidad no UCUM.")
    unit = str(unit_code.CodeValue)
    if unit not in UNITS[field]:
        raise InputError(f"{field}: unidad desconocida {unit!r}.")
    value = float(measured.NumericValue) * UNITS[field][unit]
    if not math.isfinite(value):
        raise InputError(f"{field}: número no finito.")
    return value, unit


def scalar(item: Dataset, field: str) -> tuple[float | None, list[dict]]:
    readings = []
    for node in nodes(item, CODE[field]):
        value, unit = numeric(node, field)
        readings.append({"value": value, "source_value": str(node.MeasuredValueSequence[0].NumericValue),
                         "unit": unit, "concept": CODE[field]})
    if not readings:
        return None, []
    vals = [entry["value"] for entry in readings]
    if not all(math.isclose(vals[0], v, rel_tol=1e-8, abs_tol=1e-12) for v in vals[1:]):
        raise InputError(f"{field}: medidas múltiples distintas; no se promedian silenciosamente.")
    return vals[0], readings


def textual(item: Dataset, field: str) -> str | None:
    vals = nodes(item, CODE[field])
    if not vals:
        return None
    if len(vals) > 1:
        raise InputError(f"{field}: código duplicado.")
    node = vals[0]
    vt = str(node.ValueType)
    if vt == "CODE":
        seq = node.ConceptCodeSequence[0]
        return str(seq.CodeValue)
    if vt == "UIDREF":
        return str(node.UID)
    if vt == "DATETIME":
        return str(node.DateTime)
    if vt == "TEXT":
        return str(node.TextValue)
    raise InputError(f"{field}: tipo {vt} no reconocido.")


def parse_event(item: Dataset) -> IrradiationEvent:
    orig: dict[str, list[dict]] = {}
    quality: dict[str, str] = {}
    values = {}
    for field in UNITS:
        value, readings = scalar(item, field)
        values[field] = value
        quality[field] = "present" if readings else "missing"
        if readings:
            orig[field] = readings
    uid = textual(item, "uid")
    if not uid:
        raise InputError("Irradiation Event UID (113769) ausente.")
    plane = textual(item, "plane")
    if not plane:
        raise InputError(f"Evento {uid}: Acquisition Plane (113764) ausente.")
    # Material (113757) is nested under each X-Ray Filters container (113771).
    filters = nodes(item, "113771")
    if filters:
        if len(filters) != 1:
            raise InputError(f"Evento {uid}: filtración de varios filtros; configurar modelo espectral.")
        mats = nodes(filters[0], CODE["filter_material"])
        if len(mats) != 1 or not getattr(mats[0], "ConceptCodeSequence", None):
            raise InputError(f"Evento {uid}: material filtrante indeterminado.")
        matcode = mats[0].ConceptCodeSequence[0]
        material = (str(matcode.CodingSchemeDesignator), str(matcode.CodeValue))
        if material != ("SCT", "66925006"):
            # Never reinterpret a thickness of aluminium as copper.
            raise InputError(f"Evento {uid}: filtro {material}; calibración solo admite Cu.")
    elif values["filtration_mm_cu"] is not None:
        raise InputError("Espesor de filtro sin identificación de material.")
    iso = tuple(values[k] for k in ("table_x_iso", "table_y_iso", "table_z_iso"))
    legacy = tuple(values[k] for k in ("table_lat", "table_long", "table_height"))
    if any(v is not None for v in iso) and not all(v is not None for v in iso):
        raise InputError(f"Evento {uid}: coordenadas isocéntricas incompletas.")
    if all(v is not None for v in iso):
        table, system = iso, "isocenter"
    elif all(v is not None for v in legacy):
        table, system = legacy, "legacy"
    elif any(v is not None for v in legacy):
        raise InputError(f"Evento {uid}: coordenadas de mesa incompletas.")
    else:
        table, system = None, None
    if values["exposure_mas"] is None and values["exposure_time"] is not None and values["current_ma"] is not None:
        # Exposure and average current need not multiply in pulsed operation: no substitution.
        quality["exposure_mas"] = "missing_not_inferred_from_current_time"
    elif values["exposure_mas"] is not None:
        quality["exposure_mas"] = "present"
    return IrradiationEvent(
        uid=uid, plane=plane, mode=textual(item, "mode"),
        ka_r_gy=values["ka_r_gy"], kap_gy_m2=values["kap_gy_m2"],
        kvp=values["kvp"], filtration_mm_cu=values["filtration_mm_cu"],
        exposure_mas=values["exposure_mas"], primary_deg=values["primary_deg"],
        secondary_deg=values["secondary_deg"], primary_end_deg=values["primary_end_deg"],
        secondary_end_deg=values["secondary_end_deg"], table_mm=table,
        table_system=system, sid_mm=values["sid_mm"], sod_mm=values["sod_mm"],
        srp_mm=values["srp_mm"], field_width_mm=values["field_width_mm"],
        field_height_mm=values["field_height_mm"], timestamp=textual(item, "timestamp"),
        quality=quality, original=orig)


def read_rdsr(path: str | Path) -> RDSRData:
    path = Path(path)
    dataset = pydicom.dcmread(path, stop_before_pixels=True)
    if str(getattr(dataset, "SOPClassUID", "")) != XRAY_DOSE_SR_UID:
        raise InputError("El fichero no tiene la SOP Class X-Ray Radiation Dose SR Storage.")
    event_nodes = [x for x in descendants(dataset) if concept(x) == ("DCM", EVENT)]
    if not event_nodes:
        raise InputError("RDSR sin contenedores Irradiation Event X-Ray Data (113706).")
    events = tuple(parse_event(x) for x in event_nodes)
    seen = [e.uid for e in events]
    if len(seen) != len(set(seen)):
        raise InputError("UID de evento repetido; no se duplicará su dosis.")
    totals = {}
    for k, codes, unit_field in (
        ("ka_r", ("113728", "113729"), "ka_r_gy"),
        ("kap", ("113726", "113727"), "kap_gy_m2"),
    ):
        vals = []
        for code in codes:
            for node in nodes(dataset, code):
                vals.append(numeric(node, unit_field)[0])
        totals[k] = sum(vals) if vals else None
    warnings = []
    for field, key in (("ka_r_gy", "ka_r"), ("kap_gy_m2", "kap")):
        total = totals[key]
        if total is not None:
            known = [getattr(ev, field) for ev in events]
            if all(x is not None for x in known) and not math.isclose(
                sum(known), total, abs_tol=1e-7, rel_tol=0.02
            ):
                warnings.append(f"Suma por evento de {field} difiere del total acumulado; revisar alcance/planos.")
    return RDSRData(events, str(getattr(dataset, "Manufacturer", "")),
                    str(getattr(dataset, "ManufacturerModelName", "")),
                    sha256(path.read_bytes()).hexdigest(), totals["ka_r"], totals["kap"], tuple(warnings))


def _dicom_files(root: Path) -> Iterable[Path]:
    if root.is_file():
        yield root
        return
    for path in root.rglob("*"):
        if path.is_file():
            yield path


def find_rdsr(root: str | Path) -> dict:
    """Locate X-Ray Radiation Dose SR files under a DICOM export directory.

    The scanner reads metadata only. It does not parse image pixels and does not
    export patient identifiers. XA/RF image objects are reported as non-candidates
    so a user can understand why a selected file cannot drive the calculator.
    """
    root = Path(root)
    if not root.exists():
        raise InputError(f"No existe la ruta: {root}")
    candidates = []
    rejected = Counter()
    unreadable = 0
    for path in _dicom_files(root):
        try:
            ds = pydicom.dcmread(path, stop_before_pixels=True, specific_tags=[
                "SOPClassUID", "Modality", "Manufacturer", "ManufacturerModelName",
                "StudyInstanceUID", "SeriesInstanceUID", "SOPInstanceUID"
            ])
        except Exception:
            unreadable += 1
            continue
        sop = str(getattr(ds, "SOPClassUID", ""))
        modality = str(getattr(ds, "Modality", ""))
        if sop == XRAY_DOSE_SR_UID:
            try:
                parsed = read_rdsr(path)
                event_count = len(parsed.events)
                warnings = list(parsed.warnings)
                status = "readable_rdsr"
            except Exception as exc:
                event_count = None
                warnings = [str(exc)]
                status = "rdsr_but_not_readable_by_current_parser"
            candidates.append({
                "path": str(path),
                "status": status,
                "manufacturer": str(getattr(ds, "Manufacturer", "")),
                "model": str(getattr(ds, "ManufacturerModelName", "")),
                "event_count": event_count,
                "study_uid_sha256": sha256(str(getattr(ds, "StudyInstanceUID", "")).encode()).hexdigest()
                if getattr(ds, "StudyInstanceUID", None) else None,
                "series_uid_sha256": sha256(str(getattr(ds, "SeriesInstanceUID", "")).encode()).hexdigest()
                if getattr(ds, "SeriesInstanceUID", None) else None,
                "sop_instance_uid_sha256": sha256(str(getattr(ds, "SOPInstanceUID", "")).encode()).hexdigest()
                if getattr(ds, "SOPInstanceUID", None) else None,
                "warnings": warnings,
            })
        else:
            rejected[(modality or "unknown", sop or "unknown")] += 1
    return {
        "root": str(root),
        "rdsr_count": len(candidates),
        "rdsr_candidates": candidates,
        "non_rdsr_summary": [
            {"modality": modality, "sop_class_uid": sop, "count": count}
            for (modality, sop), count in sorted(rejected.items())
        ],
        "unreadable_file_count": unreadable,
        "privacy_note": "UIDs are hashed; patient name and patient ID are not exported.",
    }


def inventory(path: str | Path) -> dict:
    """Inventory of coded names and units only; no patient identifiers or free-text values."""
    ds = pydicom.dcmread(path, stop_before_pixels=True)
    counts = Counter()
    for node in descendants(ds):
        c = concept(node)
        if c:
            units = ""
            if getattr(node, "ValueType", "") == "NUM":
                seq = getattr(node, "MeasuredValueSequence", [])
                if seq and getattr(seq[0], "MeasurementUnitsCodeSequence", []):
                    units = str(seq[0].MeasurementUnitsCodeSequence[0].CodeValue)
            counts[(c[0], c[1], str(node.ValueType), units)] += 1
    return {"sop_class": str(getattr(ds, "SOPClassUID", "")), "concepts": [
        {"scheme": a, "code": b, "value_type": c, "unit": d, "count": n}
        for (a, b, c, d), n in sorted(counts.items())]}
