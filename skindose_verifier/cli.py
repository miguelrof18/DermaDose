"""Interfaz reproducible: inspección, cálculo, simulación y validación."""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import replace
from hashlib import sha256
from importlib.resources import files
import json
from pathlib import Path
import shutil
import sys

from . import __version__
from .config import Equipment
from .convergence import compare_resolutions
from .domain import InputError
from .engine import calculate
from .rdsr import find_rdsr, inventory, read_rdsr
from .reporting import save_results
from .sensitivity import position_sensitivity
from .synthetic import make_synthetic_rdsr
from .validation import compare_with_radimetrics, validate_measurements, validate_map
from .provenance import runtime_provenance
from .uncertainty import propagate_parameters
from .benchmark import benchmark_calculation
from .mesh_processing import mesh_quality
from .patient import sample_patient
import numpy as np
from scipy.spatial import cKDTree


def run(rdsr: Path, config: Path, stage: str, destination: Path) -> dict:
    equipment = Equipment.load(config)
    report = read_rdsr(rdsr)
    events = report.events
    synthetic_override = False
    try:
        equipment.require_matching_device(report.manufacturer, report.model)
    except InputError:
        if equipment.raw.get("azurion_profile") or equipment.raw.get("allow_synthetic_device_override") is not True:
            raise
        synthetic_override = True
        if equipment.surface_factor_provider is not None:
            equipment.surface_factor_provider.spec["field_shape_policy"] = "equal_area"
        # One configured apparatus profile: retain event kerma, angles and beam
        # quality; use configured geometry where the foreign RDSR frame cannot
        # be mapped. Record every assumption in the event and global audit.
        converted = []
        for event in events:
            geometry = equipment.geometry_for(event.plane)
            quality = {**event.quality,
                       "table_mm": "assumed_fixed_from_configuration",
                       "sid_mm": "assumed_from_configuration",
                       "sod_mm": "assumed_from_configuration",
                       "srp_mm": "assumed_from_configuration",
                       "synthetic_profile_spectrum": "additional_al_not_modeled"}
            field_width = event.field_width_mm
            field_height = event.field_height_mm
            shutter = event.original.get("philips_shutters_mm_at_1000_mm", {})
            if field_width is None and field_height is None and set(shutter) == {"left", "right", "top", "bottom"}:
                scale = float(geometry["source_to_detector_mm"]) / 1000.0
                field_width = (shutter["left"] + shutter["right"]) * scale
                field_height = (shutter["top"] + shutter["bottom"]) * scale
                quality["field_size"] = "philips_shutters_scaled_from_1000_mm_to_configured_detector"
            converted.append(replace(event, table_mm=tuple(geometry["table_origin_mm"]),
                table_system=geometry["table_coordinates"],
                sid_mm=float(geometry["source_to_detector_mm"]),
                sod_mm=float(geometry["source_to_isocenter_mm"]),
                srp_mm=float(geometry["source_to_irp_mm"]),
                field_width_mm=field_width, field_height_mm=field_height, quality=quality))
        events = tuple(converted)
    factor_source = "configured"
    fallback_reason = None
    try:
        calculation = calculate(events, equipment, stage)
    except InputError as exc:
        if equipment.raw.get("azurion_profile") or equipment.raw.get("allow_historical_factor_fallback") is not True:
            raise
        if (stage != "skin_dose" or equipment.raw.get("surface_factors", {}).get("provider") != "benmakhlouf_2011"
                or "fuera del dominio" not in str(exc)):
            raise
        # A whole-run switch avoids mixing published and synthetic factors.
        # The source RDSR and the user's equipment JSON remain unchanged.
        synthetic = json.loads(files("skindose_verifier").joinpath(
            "data/equipo_sintetico.json").read_text(encoding="utf-8"))["factors"]
        replacement = deepcopy(equipment.raw)
        replacement.pop("surface_factors", None)
        replacement["factors"]["bsf"] = synthetic["bsf"]
        replacement["factors"]["f_skin"] = synthetic["f_skin"]
        equipment = Equipment.from_dict(replacement)
        calculation = calculate(events, equipment, stage)
        factor_source = "historical_synthetic_whole_run"
        fallback_reason = str(exc)
        calculation.summary["status"] = "research_synthetic_factors_out_of_domain_no_clinical_use"
    calculation.summary["factor_source"] = factor_source
    audit = {"software_version": __version__, "rdsr_sha256": report.sha256,
             "config_sha256": sha256(config.read_bytes()).hexdigest(),
             "mesh_sha256": sha256(Path(equipment.patient["mesh_path"]).read_bytes()).hexdigest()
             if equipment.patient["type"] == "mesh" else None,
             "configuration_version": equipment.version,
             "azurion_profile": equipment.raw.get("azurion_profile"),
             "validation_status": equipment.raw.get("validation_status", "not_reported"),
             "synthetic_equipment_profile_used_for_different_rdsr_device": synthetic_override,
             "effective_factor_source": factor_source,
             "factor_fallback_reason": fallback_reason,
             "rdsr_device": {"manufacturer": report.manufacturer, "model": report.model},
             "configuration_device": {"manufacturer": equipment.raw["manufacturer"],
                                      "model": equipment.raw["model"]},
             "phantom": equipment.patient.get("phantom_name", "XCAT50"),
             "accumulated_ka_r_gy": report.accumulated_ka_r_gy,
             "accumulated_kap_gy_m2": report.accumulated_kap_gy_m2,
             "rdsr_warnings": report.warnings, "runtime": runtime_provenance()}
    return {"files": save_results(calculation, destination, audit=audit),
            "quantity": calculation.summary["quantity"],
            "psd_point_gy": calculation.summary["psd_point_gy"],
            "status": calculation.summary["status"]}


def compare_phantoms(rdsr: Path, xcat_config: Path, alternative_config: Path,
                     stage: str, destination: Path) -> dict:
    """Run one RDSR on two phantom profiles and write an anatomy comparison."""
    destination.mkdir(parents=True, exist_ok=True)
    xcat = run(rdsr, xcat_config, stage, destination / "XCAT50")
    alt = run(rdsr, alternative_config, stage, destination / "alternative")
    xj = json.loads(Path(xcat["files"]["summary"]).read_text(encoding="utf-8"))
    aj = json.loads(Path(alt["files"]["summary"]).read_text(encoding="utf-8"))
    xpsd = float(xj["psd_point_gy"])
    apsd = float(aj["psd_point_gy"])
    relative = None if xpsd == 0 else (apsd - xpsd) / xpsd * 100.0
    def metrics(summary):
        m = summary["surface_metrics"]["irradiated_surface"]
        return {"psd_gy": summary["psd_point_gy"], "dmean_gy": m["dmean_gy"],
                "dmax_gy": m["dmax_gy"], "d2_percent_gy": m["d2_percent_gy"],
                "d5_percent_gy": m["d5_percent_gy"], "d10_percent_gy": m["d10_percent_gy"],
                "irradiated_area_cm2": m["area_cm2"],
                "dose_area_integral_gy_cm2": summary["surface_metrics"]["dose_area_integral_gy_cm2"],
                "event_count": summary["event_count"]}
    comparison = {"baseline": "XCAT50", "alternative": json.loads(alternative_config.read_text(encoding="utf-8")).get("phantom", {}).get("name", "alternative"),
                  "stage": stage, "xcat50": metrics(xj), "alternative": metrics(aj),
                  "relative_psd_difference_percent": relative,
                  "algorithm_constant": "same RDSR, factors, geometry and interpolation; only patient mesh/profile changes"}
    xn = np.load(xcat["files"]["numeric_map"]); an = np.load(alt["files"]["numeric_map"])
    xp, ap = xn["patient_points_mm"], an["patient_points_mm"]
    xd, ad = xn["dose_gy"].ravel(), an["dose_gy"].ravel()
    tree = cKDTree(xp); dist, idx = tree.query(ap)
    comparison["spatial_difference_on_alternative_points"] = {
        "nearest_xcat_distance_mm_mean": float(np.mean(dist)),
        "nearest_xcat_distance_mm_max": float(np.max(dist)),
        "dose_difference_gy_mean": float(np.mean(ad - xd[idx])),
        "dose_difference_gy_rms": float(np.sqrt(np.mean((ad - xd[idx]) ** 2))),
        "note": "nearest-surface comparison; meshes are not pointwise identical"}
    (destination / "comparacion_fantomas.json").write_text(
        json.dumps(comparison, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    return comparison


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="skindose", description="Prototipo de investigación, NO de uso clínico")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    i = sub.add_parser("inspect", help="Inventario de códigos/unidades DICOM sin identificadores personales")
    i.add_argument("rdsr", type=Path)
    fr = sub.add_parser("find-rdsr", help="Buscar RDSR dentro de una carpeta DICOM o comprobar un fichero")
    fr.add_argument("path", type=Path)
    comp = sub.add_parser("compatibility-report", help="Informe de compatibilidad de RDSR reales anonimizados")
    comp.add_argument("path", type=Path, help="RDSR o carpeta DICOM")
    comp.add_argument("--config", type=Path, help="Configuración de equipo para comprobar fabricante/modelo")
    comp.add_argument("--stage", choices=["free_air", "calibrated_air", "skin_dose"], default="skin_dose")
    d = sub.add_parser("demo", help="Generar RDSR totalmente sintético y mapa de ejemplo")
    d.add_argument("--out", type=Path, default=Path("demo_output"))
    d.add_argument("--model", choices=["plane", "ellipse", "mesh", "mesh-pyskindose"], default="plane")
    d.add_argument("--stage", choices=["free_air", "calibrated_air", "skin_dose"], default="free_air")
    d.add_argument("--mesh-edge-mm", type=float, help="Subdividir aristas conservando el XCAT original")
    d.add_argument("--mesh-roi", help="Refinamiento regional: xmin,ymin,zmin,xmax,ymax,zmax, en mm de paciente")
    b = sub.add_parser("run", help="Calcular sobre un RDSR con configuración explícita")
    b.add_argument("rdsr", type=Path)
    b.add_argument("--config", type=Path, required=True)
    b.add_argument("--out", type=Path, required=True)
    b.add_argument("--stage", choices=["free_air", "calibrated_air", "skin_dose"], default="free_air")
    pc = sub.add_parser("compare-phantoms", help="Ejecutar el mismo RDSR en XCAT50 y otro fantoma")
    pc.add_argument("rdsr", type=Path)
    pc.add_argument("--xcat-config", type=Path, required=True)
    pc.add_argument("--alternative-config", type=Path, required=True)
    pc.add_argument("--out", type=Path, required=True)
    pc.add_argument("--stage", choices=["free_air", "calibrated_air", "skin_dose"], default="skin_dose")
    q = sub.add_parser("convergence", help="Comprobar sensibilidad de PSD al paso de la malla")
    q.add_argument("rdsr", type=Path)
    q.add_argument("--config", type=Path, required=True)
    q.add_argument("--resolutions", default="10,5,2.5", help="Resoluciones en mm, separadas por comas")
    q.add_argument("--stage", choices=["free_air", "calibrated_air", "skin_dose"], default="skin_dose")
    q.add_argument("--out", type=Path)
    u = sub.add_parser("uncertainty", help="Propagación de parámetros LHS/Monte Carlo; no transporte de fotones")
    u.add_argument("rdsr", type=Path)
    u.add_argument("--config", type=Path, required=True)
    u.add_argument("--spec", type=Path, required=True)
    u.add_argument("--out", type=Path, required=True, help="Directorio JSON/NPZ")
    u.add_argument("--stage", choices=["free_air", "calibrated_air", "skin_dose"], default="skin_dose")
    bench = sub.add_parser("benchmark", help="Comparar cálculo con y sin caché geométrica")
    bench.add_argument("rdsr", type=Path)
    bench.add_argument("--config", type=Path, required=True)
    bench.add_argument("--out", type=Path, required=True)
    bench.add_argument("--repeats", type=int, default=3)
    bench.add_argument("--fields", type=int, default=12)
    bench.add_argument("--stage", choices=["free_air", "calibrated_air", "skin_dose"], default="skin_dose")
    mq = sub.add_parser("mesh-quality", help="Auditar geometría original o refinada")
    mq.add_argument("--config", type=Path, required=True)
    mq.add_argument("--out", type=Path, required=True)
    s = sub.add_parser("sensitivity-position", help="Recalcular con posición de paciente desplazada ±mm por eje")
    s.add_argument("rdsr", type=Path)
    s.add_argument("--config", type=Path, required=True)
    s.add_argument("--offset-mm", type=float, default=5.)
    s.add_argument("--stage", choices=["free_air", "calibrated_air", "skin_dose"], default="skin_dose")
    s.add_argument("--out", type=Path, help="JSON opcional para conservar el estudio de sensibilidad")
    v = sub.add_parser("validate-measurements", help="Comparar con medidas independientes de fantoma")
    v.add_argument("cases", type=Path)
    v.add_argument("--stage", choices=["free_air", "calibrated_air", "skin_dose"], required=True)
    m = sub.add_parser("validate-map", help="Comparar mapa con puntos medidos ya registrados")
    m.add_argument("results", type=Path, help="resultados.json")
    m.add_argument("measurements", type=Path, help="CSV de puntos experimentales")
    m.add_argument("--threshold", type=float, default=.5, help="Fracción del máximo medido para isodosis")
    m.add_argument("--max-distance-mm", type=float, default=5.)
    c = sub.add_parser("compare-radimetrics", help="Comparación secundaria, sin validar físicamente el modelo")
    c.add_argument("results", type=Path, help="resultados.json del calculador")
    c.add_argument("--kar", type=float, required=True)
    c.add_argument("--psd", type=float, required=True)
    c.add_argument("--medium", choices=["skin", "water"], required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "inspect":
            output = inventory(args.rdsr)
        elif args.command == "find-rdsr":
            output = find_rdsr(args.path)
        elif args.command == "compatibility-report":
            output = compatibility_report(args.path, args.config, args.stage)
        elif args.command == "demo":
            args.out.mkdir(parents=True, exist_ok=True)
            config = args.out/"equipo_sintetico.json"
            template = files("skindose_verifier").joinpath("data/equipo_sintetico.json")
            data = json.loads(template.read_text(encoding="utf-8"))
            if args.model == "ellipse":
                data["patient"] = {"type": "ellipse", "semi_axis_x_mm": 170,
                                   "semi_axis_z_mm": 110, "length_mm": 500,
                                   "origin_on_table_mm": [0, 0, 0], "resolution_mm": 5}
            elif args.model == "mesh":
                asset = "xcat_wholebody_body.stl"
                shutil.copyfile(files("skindose_verifier").joinpath("data", asset), args.out/asset)
                data["patient"] = {"type": "mesh", "mesh_path": asset, "mesh_unit": "mm",
                                   "mesh_axes_matrix": [[1, 0, 0], [0, 0, 1], [0, -1, 0]],
                                   "mesh_origin_mm": [0, -500, 40], "thickness_mm": 220,
                                   "origin_on_table_mm": [0, 0, 0], "resolution_mm": 5}
            elif args.model == "mesh-pyskindose":
                asset = "pyskindose_adult_female.stl"
                shutil.copyfile(files("skindose_verifier").joinpath("data", asset), args.out/asset)
                data["patient"] = {"type": "mesh", "mesh_path": asset, "mesh_unit": "cm",
                                   "mesh_axes_matrix": [[1, 0, 0], [0, 0, 1], [0, -1, 0]],
                                   "mesh_origin_mm": [0, 350, -100], "thickness_mm": 220,
                                   "origin_on_table_mm": [0, 0, 0], "resolution_mm": 5}
            if args.mesh_edge_mm is not None:
                if data["patient"]["type"] != "mesh":
                    raise InputError("--mesh-edge-mm requiere modelo mesh.")
                data["patient"]["mesh_max_edge_mm"] = args.mesh_edge_mm
            if args.mesh_roi is not None:
                if args.mesh_edge_mm is None:
                    raise InputError("--mesh-roi requiere --mesh-edge-mm.")
                bounds = [float(x) for x in args.mesh_roi.split(",")]
                if len(bounds) != 6:
                    raise InputError("--mesh-roi requiere seis coordenadas.")
                data["patient"]["mesh_refinement_roi_mm"] = [bounds[:3], bounds[3:]]
            config.write_text(json.dumps(data, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
            rdsr = make_synthetic_rdsr(args.out/"rdsr_sintetico.dcm")
            output = run(rdsr, config, args.stage, args.out)
        elif args.command == "run":
            output = run(args.rdsr, args.config, args.stage, args.out)
        elif args.command == "compare-phantoms":
            output = compare_phantoms(args.rdsr, args.xcat_config, args.alternative_config,
                                      args.stage, args.out)
        elif args.command == "convergence":
            equipment = Equipment.load(args.config)
            report = read_rdsr(args.rdsr)
            equipment.require_matching_device(report.manufacturer, report.model)
            output = compare_resolutions(report.events, equipment, args.stage,
                                         [float(x.strip()) for x in args.resolutions.split(",")])
            output["runtime"] = runtime_provenance()
            if args.out:
                _write_json(args.out, output)
        elif args.command in ("uncertainty", "benchmark"):
            equipment = Equipment.load(args.config)
            report = read_rdsr(args.rdsr)
            equipment.require_matching_device(report.manufacturer, report.model)
            if args.command == "benchmark":
                output = benchmark_calculation(report.events, equipment, args.stage, args.repeats, args.fields)
                _write_json(args.out, output)
            else:
                spec = json.loads(args.spec.read_text(encoding="utf-8"))
                propagated = propagate_parameters(report.events, equipment, args.stage, spec)
                output = propagated.summary
                output["audit"] = {"rdsr_sha256": report.sha256, "config_sha256": sha256(args.config.read_bytes()).hexdigest(),
                                   "runtime": runtime_provenance()}
                _write_json(args.out/"incertidumbre.json", output)
                np.savez_compressed(args.out/"mapa_incertidumbre.npz", patient_points_mm=propagated.points_mm,
                                    area_mm2=propagated.area_mm2, mean_gy=propagated.mean_gy,
                                    standard_deviation_gy=propagated.standard_deviation_gy,
                                    probability_exceeding_threshold=propagated.probability_exceeding_threshold,
                                    threshold_gy=output["threshold_gy"])
        elif args.command == "mesh-quality":
            surface = sample_patient(Equipment.load(args.config).patient)
            if surface.mesh is None:
                raise InputError("mesh-quality requiere un modelo mesh.")
            output = mesh_quality(surface.mesh)
            _write_json(args.out, output)
        elif args.command == "sensitivity-position":
            equipment = Equipment.load(args.config)
            report = read_rdsr(args.rdsr)
            equipment.require_matching_device(report.manufacturer, report.model)
            output = position_sensitivity(report.events, equipment, args.stage, args.offset_mm)
            output["audit"] = {"software_version": __version__, "rdsr_sha256": report.sha256,
                               "config_sha256": sha256(args.config.read_bytes()).hexdigest()}
            if args.out is not None:
                args.out.parent.mkdir(parents=True, exist_ok=True)
                args.out.write_text(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False)+"\n",
                                    encoding="utf-8")
        elif args.command == "validate-measurements":
            output = validate_measurements(args.cases, args.stage)
        elif args.command == "validate-map":
            output = validate_map(args.results, args.measurements, args.threshold, args.max_distance_mm)
        else:
            ours = json.loads(args.results.read_text(encoding="utf-8"))
            output = compare_with_radimetrics(ours, args.kar, args.psd, args.medium)
    except (InputError, OSError, ValueError, KeyError) as exc:
        print(f"ERROR CONTROLADO: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


def compatibility_report(path: Path, config: Path | None, stage: str) -> dict:
    found = find_rdsr(path)
    equipment = Equipment.load(config) if config else None
    reports = []
    required = ("ka_r_gy", "primary_deg", "secondary_deg", "table_mm", "sid_mm",
                "sod_mm", "srp_mm", "field_width_mm", "field_height_mm")
    for candidate in found["rdsr_candidates"]:
        item = {**candidate}
        try:
            report = read_rdsr(candidate["path"])
            item["readable_by_parser"] = True
            item["manufacturer"] = report.manufacturer
            item["model"] = report.model
            item["event_count"] = len(report.events)
            item["event_missing_fields"] = {
                ev.uid: [field for field in required if getattr(ev, field) is None]
                for ev in report.events
            }
            item["events_with_missing_required_fields"] = sum(
                1 for missing in item["event_missing_fields"].values() if missing
            )
            if stage != "free_air":
                item["events_missing_beam_quality"] = sum(
                    1 for ev in report.events if ev.kvp is None or ev.filtration_mm_cu is None
                )
            if equipment is not None:
                try:
                    equipment.require_matching_device(report.manufacturer, report.model)
                    item["configuration_device_match"] = True
                    item["configuration_device_message"] = "manufacturer_and_model_match"
                except InputError as exc:
                    item["configuration_device_match"] = False
                    item["configuration_device_message"] = str(exc)
            item["compatible_for_stage_without_full_dose_run"] = (
                item["events_with_missing_required_fields"] == 0 and
                (stage == "free_air" or item.get("events_missing_beam_quality", 0) == 0) and
                (equipment is None or item.get("configuration_device_match") is True)
            )
        except Exception as exc:
            item["readable_by_parser"] = False
            item["compatible_for_stage_without_full_dose_run"] = False
            item["error"] = str(exc)
        reports.append(item)
    return {
        "input": str(path),
        "stage": stage,
        "rdsr_count": found["rdsr_count"],
        "non_rdsr_summary": found["non_rdsr_summary"],
        "unreadable_file_count": found["unreadable_file_count"],
        "compatibility_reports": reports,
        "privacy_note": "PatientName and PatientID are not exported; UIDs are hashed in discovery output.",
    }


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
