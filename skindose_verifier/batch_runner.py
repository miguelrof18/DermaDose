"""Run DermaDose on one folder per patient, without interactive prompts."""
import argparse
import json
from pathlib import Path
import subprocess
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    config = args.config.resolve()
    if not root.is_dir() or not config.is_file():
        parser.error("Indique un directorio de casos y un JSON de equipo existentes.")
    cases = sorted(p for p in root.iterdir() if p.is_dir())
    summary = []
    for case in cases:
        command = [sys.executable, "-m", "skindose_verifier.example_runner", str(case),
                   "--mesh", "patient", "--config", str(config)]
        if args.prepare_only:
            command.append("--prepare-only")
        audit_file = case / "resultados" / "seleccion_entrada.json"
        previous = audit_file.stat().st_mtime_ns if audit_file.exists() else None
        print(f"Procesando {case.name}", flush=True)
        with (case / "dermadose_lote.log").open("w", encoding="utf-8") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
        audit = {}
        if audit_file.exists() and audit_file.stat().st_mtime_ns != previous:
            audit = json.loads(audit_file.read_text(encoding="utf-8"))
        summary.append(dict(case=case.name, exit_code=result.returncode,
                            status=("preparado" if args.prepare_only else "calculado")
                            if result.returncode == 0 else
                            ("informacion_insuficiente" if audit.get("registration_status")
                             == "insufficient_information" else "error"),
                            placement=audit.get("registration_status"),
                            reason=audit.get("registration_note"),
                            log=str(case / "dermadose_lote.log")))
        (root / "dermadose_lote.json").write_text(
            json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return 0 if cases and all(row["exit_code"] == 0 for row in summary) else 1


if __name__ == "__main__":
    raise SystemExit(main())
