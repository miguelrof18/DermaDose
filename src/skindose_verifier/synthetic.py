"""RDSR artificial para pruebas. No deriva de un paciente ni un equipo clínico."""
from pathlib import Path

from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import ExplicitVRLittleEndian

from .rdsr import XRAY_DOSE_SR_UID

ROOT_UID = "1.2.826.0.1.3680043.10.543.99"


def coded(code: str, meaning: str, scheme: str = "DCM") -> Dataset:
    ds = Dataset()
    ds.CodeValue, ds.CodingSchemeDesignator, ds.CodeMeaning = code, scheme, meaning
    return ds


def item(value_type: str, code: str, name: str, *, value=None, unit: str | None = None,
         value_code: Dataset | None = None) -> Dataset:
    ds = Dataset()
    ds.ValueType = value_type
    ds.RelationshipType = "CONTAINS"
    ds.ConceptNameCodeSequence = Sequence([coded(code, name)])
    if value_type == "CONTAINER":
        ds.ContinuityOfContent = "SEPARATE"
        ds.ContentSequence = Sequence([])
    elif value_type == "NUM":
        m = Dataset()
        m.NumericValue = str(value)
        m.MeasurementUnitsCodeSequence = Sequence([coded(unit, unit, "UCUM")])
        ds.MeasuredValueSequence = Sequence([m])
    elif value_type == "CODE":
        ds.ConceptCodeSequence = Sequence([value_code])
    elif value_type == "UIDREF":
        ds.UID = str(value)
    elif value_type == "DATETIME":
        ds.DateTime = str(value)
    return ds


def make_synthetic_rdsr(path: str | Path, *, two_events: bool = True,
                        omit: str | None = None, biplane: bool = False) -> Path:
    """Fully fictional fixture; a self-consistency check, never clinical validation."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fm = FileMetaDataset()
    fm.MediaStorageSOPClassUID = XRAY_DOSE_SR_UID
    fm.MediaStorageSOPInstanceUID = ROOT_UID + ".1"
    fm.TransferSyntaxUID = ExplicitVRLittleEndian
    fm.ImplementationClassUID = ROOT_UID + ".2"
    ds = FileDataset(str(path), {}, file_meta=fm, preamble=b"\x00"*128)
    ds.SOPClassUID, ds.SOPInstanceUID = fm.MediaStorageSOPClassUID, fm.MediaStorageSOPInstanceUID
    ds.Modality = "SR"
    ds.Manufacturer, ds.ManufacturerModelName = "SIMULATED", "PHANTOM-ONLY"
    ds.PatientName = "SYNTHETIC^NO_PERSON"
    ds.PatientID = "SYNTHETIC_ONLY"
    ds.ContentDate, ds.ContentTime = "20250101", "120000"
    root = item("CONTAINER", "113701", "X-Ray Radiation Dose Report")
    del root.RelationshipType
    ds.ContentSequence = Sequence([root])
    n = 2 if two_events else 1
    for idx in range(n):
        ev = item("CONTAINER", "113706", "Irradiation Event X-Ray Data")
        root.ContentSequence.append(ev)
        plane = "113621" if (idx == 1 and biplane) else "113620" if biplane else "113622"
        ev.ContentSequence.extend([
            item("CODE", "113764", "Acquisition Plane", value_code=coded(plane, "Plane")),
            item("UIDREF", "113769", "Irradiation Event UID", value=ROOT_UID+f".10.{idx+1}"),
            item("DATETIME", "111526", "DateTime Started", value=f"2025010112000{idx}"),
            item("CODE", "113721", "Irradiation Event Type",
                 value_code=coded("44491008", "Fluoroscopy", "SCT")),
        ])
        numbers = [
            ("113738", "Dose (RP)", .019834710744, "Gy"),
            ("122130", "Dose Area Product", .00019834710744, "Gy.m2"),
            ("113733", "KVP", 80, "kV"),
            ("113736", "Exposure", 20000, "uA.s"),
            ("112011", "Positioner Primary Angle", 0 if idx == 0 else 10, "deg"),
            ("112012", "Positioner Secondary Angle", 0, "deg"),
            ("113750", "Distance Source to Detector", 1100, "mm"),
            ("113748", "Distance Source to Isocenter", 700, "mm"),
            ("113737", "Distance Source to Reference Point", 550, "mm"),
            ("113789", "Collimated Field Width", 200, "mm"),
            ("113788", "Collimated Field Height", 200, "mm"),
            ("128766", "Table X Position to Isocenter", 0, "mm"),
            ("128767", "Table Y Position to Isocenter", 0, "mm"),
            ("128768", "Table Z Position to Isocenter", 0, "mm"),
        ]
        for code, name, value, unit in numbers:
            if code != omit:
                ev.ContentSequence.append(item("NUM", code, name, value=value, unit=unit))
        filters = item("CONTAINER", "113771", "X-Ray Filters")
        filters.ContentSequence.extend([
            item("CODE", "113757", "X-Ray Filter Material",
                 value_code=coded("66925006", "Copper", "SCT")),
            item("NUM", "113758", "X-Ray Filter Thickness Minimum", value=.2, unit="mm")])
        ev.ContentSequence.append(filters)
    root.ContentSequence.extend([
        item("NUM", "113728", "Fluoro Dose (RP) Total", value=.019834710744*n, unit="Gy"),
        item("NUM", "113726", "Fluoro Dose Area Product Total", value=.00019834710744*n, unit="Gy.m2")])
    ds.save_as(path, enforce_file_format=True)
    return path
