"""Compare point PSD on two independently simplified meshes.

Usage: python -m skindose_verifier.mesh_convergence RDSR --config equipment.json
       --low surface_150k.stl --high surface_300k.stl --out comparison_psd.json
"""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import numpy as np
from .config import Equipment
from .engine import calculate
from .rdsr import read_rdsr

def compare(rdsr, config, low, high, stage="skin_dose"):
    config = Path(config).resolve()
    raw = json.loads(config.read_text(encoding="utf-8"))
    if raw.get("patient", {}).get("type") != "mesh" or raw.get("phantom"):
        raise ValueError("Utilice una configuración de malla de paciente, con su registro ya revisado.")
    report = read_rdsr(Path(rdsr))
    output = {}
    for label, path in (("150k", low), ("300k", high)):
        path = Path(path).resolve()
        data = deepcopy(raw)
        data["patient"]["mesh_path"] = str(path)
        equipment = Equipment.from_dict(data, base_directory=config.parent)
        equipment.require_matching_device(report.manufacturer, report.model)
        calculation = calculate(report.events, equipment, stage, detailed=False)
        output[label] = dict(psd_point_gy=calculation.summary["psd_point_gy"],
            peak_patient_mm=calculation.summary["psd_point_patient_mm"],
            mesh_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    a, b = output["150k"], output["300k"]
    delta = a["psd_point_gy"] - b["psd_point_gy"]
    return dict(meshes=output, stage=stage, difference_gy=delta,
        difference_percent_of_300k=None if b["psd_point_gy"] == 0 else delta/b["psd_point_gy"]*100,
        peak_displacement_mm=float(np.linalg.norm(
            np.asarray(a["peak_patient_mm"])-np.asarray(b["peak_patient_mm"]))),
        rdsr_sha256=hashlib.sha256(Path(rdsr).read_bytes()).hexdigest(),
        config_sha256=hashlib.sha256(config.read_bytes()).hexdigest(),
        precision_validated=False,
        interpretation="Una comparación de convergencia; 300k no es una solución exacta. "
                       "Verifique que ambas mallas proceden del mismo original y registro.")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rdsr",type=Path)
    parser.add_argument("--config",type=Path,required=True)
    parser.add_argument("--low",type=Path,required=True)
    parser.add_argument("--high",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--stage",choices=["free_air","calibrated_air","skin_dose"],default="skin_dose")
    args=parser.parse_args()
    if args.out.exists():
        parser.error("La salida ya existe.")
    result=compare(args.rdsr,args.config,args.low,args.high,args.stage)
    with args.out.open("x",encoding="utf-8") as handle:
        json.dump(result,handle,indent=2,ensure_ascii=False,allow_nan=False)
        handle.write("\n")
    print(json.dumps(result,indent=2,ensure_ascii=False))
if __name__ == "__main__":
    main()
