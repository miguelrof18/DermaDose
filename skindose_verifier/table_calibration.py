"""Table motion matrix S from single-axis moves recorded in an RDSR.

Implements S_.j = Delta q_I / Delta m_j (docs/sistemas_de_referencia.md, section 03)
for the shots 4-6 of the «Calibración geométrica» sheet: from a frontal position,
move only one table axis, note where the board moved (seen from the supine
patient) and, if measured with a fixed room mark, how far.

    python -m skindose_verifier.table_calibration serie.dcm \\
        --perfil config/mesa_azurion_hemodinamica_1.json \\
        --movimiento 4:izquierda:100 --movimiento 5:arriba --movimiento 6:cabecera:150 \\
        --out mesa_medida.json

EVENT is the 1-based event number in the RDSR. Each move is compared with the
previous event unless written BASE>EVENT (e.g. 3>4). Without a measured travel,
unit gain (|Delta q| = |Delta m|) is assumed and recorded as such.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path

import numpy as np

from .domain import InputError
from .placement import rdsr_readings
from .reference_frames import READING_NAMES, _matrix_with_unknown

# Directions of the board in I (x left, y down, z head) for a supine, head-first patient.
DIRECTIONS = {"izquierda": (1., 0., 0.), "derecha": (-1., 0., 0.), "abajo": (0., 1., 0.),
              "arriba": (0., -1., 0.), "cabecera": (0., 0., 1.), "pies": (0., 0., -1.),
              "left": (1., 0., 0.), "right": (-1., 0., 0.), "down": (0., 1., 0.), "up": (0., -1., 0.),
              "head": (0., 0., 1.), "feet": (0., 0., -1.)}
UNCHANGED_TOLERANCE_MM = 0.05


def parse_move(text: str) -> dict:
    """'[BASE>]EVENT:DIRECTION[:TRAVEL_MM]'."""
    parts = text.split(":")
    if len(parts) not in (2, 3):
        raise InputError(f"Movimiento {text!r}: use EVENTO:SENTIDO[:RECORRIDO_MM].")
    base, _, event = parts[0].rpartition(">")
    try:
        move = {"event": int(event), "baseline": int(base) if base else None,
                "direction": parts[1].strip().casefold(),
                "travel_mm": float(parts[2].replace(",", ".")) if len(parts) == 3 else None}
    except ValueError as exc:
        raise InputError(f"Movimiento {text!r}: números de evento o recorrido inválidos.") from exc
    if move["direction"] not in DIRECTIONS:
        raise InputError(f"Sentido {parts[1]!r} desconocido; use {', '.join(list(DIRECTIONS)[:6])}.")
    return move


def measure_columns(events, moves) -> tuple[dict[int, np.ndarray], list[dict]]:
    """One column of S per move; exactly one reading must change in each move."""
    columns, audit = {}, []
    for move in moves:
        k = move["event"]
        base = move["baseline"] if move["baseline"] is not None else k - 1
        if not (1 <= base <= len(events) and 1 <= k <= len(events)) or base == k:
            raise InputError(f"Movimiento al evento {k}: eventos fuera del RDSR ({len(events)} eventos).")
        before, after = rdsr_readings(events[base - 1]), rdsr_readings(events[k - 1])
        if before is None or after is None:
            raise InputError(f"Eventos {base} y {k}: faltan lecturas DCM 113751–113753.")
        delta = after - before
        moved = np.flatnonzero(np.abs(delta) > UNCHANGED_TOLERANCE_MM)
        if len(moved) != 1:
            raise InputError(f"Eventos {base}→{k}: cambian {len(moved)} lecturas ({delta.round(2).tolist()} mm); "
                             "mueva un solo eje de la mesa en cada toma.")
        j = int(moved[0])
        if j in columns:
            raise InputError(f"La lectura {READING_NAMES[j]} se mide dos veces.")
        travel = abs(delta[j]) if move["travel_mm"] is None else float(move["travel_mm"])
        if not np.isfinite(travel) or travel <= 0:
            raise InputError("El recorrido medido debe ser positivo.")
        physical = np.asarray(DIRECTIONS[move["direction"]]) * travel
        columns[j] = physical / delta[j]
        audit.append({"reading": READING_NAMES[j], "events": [base, k],
                      "event_uids": [events[base - 1].uid, events[k - 1].uid],
                      "delta_long_lat_height_mm": delta.tolist(), "observed_direction": move["direction"],
                      "travel_mm": travel, "travel_source": "measured" if move["travel_mm"] is not None
                      else "unit_gain_assumed", "column_i_per_mm": columns[j].tolist()})
    return columns, audit


def calibrated_profile(profile: dict, events, moves) -> dict:
    """Copy of a table profile with the measured columns of S replaced."""
    columns, audit = measure_columns(events, moves)
    result = deepcopy(profile)
    matrix = _matrix_with_unknown(profile.get("motion_matrix_i_per_mm", [[None] * 3] * 3),
                                  "S (motion_matrix_i_per_mm)")
    status = dict(profile.get("motion_matrix_status") or {})
    for j, column in columns.items():
        matrix[:, j] = column
        status[READING_NAMES[j]] = ("measured_rdsr_moves" if audit[[a["reading"] for a in audit].index(READING_NAMES[j])]
                                    ["travel_source"] == "measured" else "measured_direction_unit_gain")
    result.pop("axis_trials", None)
    result["motion_matrix_i_per_mm"] = [[None if not np.isfinite(v) else round(float(v), 9) for v in row] for row in matrix]
    result["motion_matrix_status"] = status
    result["motion_matrix_evidence"] = audit
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("rdsr", type=Path)
    parser.add_argument("--perfil", type=Path, required=True, help="Perfil JSON de la mesa de partida")
    parser.add_argument("--movimiento", action="append", required=True, help="[BASE>]EVENTO:SENTIDO[:RECORRIDO_MM]")
    parser.add_argument("--out", type=Path, help="Nuevo perfil de mesa; no se sobrescribe un archivo existente")
    args = parser.parse_args(argv)
    try:
        from .rdsr import read_rdsr
        events = read_rdsr(args.rdsr).events
        profile = json.loads(args.perfil.read_text(encoding="utf-8"))
        result = calibrated_profile(profile, events, [parse_move(m) for m in args.movimiento])
        text = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
        if args.out:
            with args.out.open("x", encoding="utf-8") as handle:
                handle.write(text)
        print(json.dumps({"motion_matrix_i_per_mm": result["motion_matrix_i_per_mm"],
                          "motion_matrix_status": result["motion_matrix_status"]}, indent=2, ensure_ascii=False))
    except (InputError, OSError, ValueError) as exc:
        parser.exit(2, f"ERROR CONTROLADO: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
