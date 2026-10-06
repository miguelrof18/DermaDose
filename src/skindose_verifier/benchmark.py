"""Comparable CPU benchmarks; export time is deliberately excluded."""
from dataclasses import replace
import time
import tracemalloc

import numpy as np

from .domain import InputError
from .engine import calculate
from .patient import sample_patient
from .provenance import runtime_provenance


def benchmark_calculation(events, equipment, stage, repeats=3, repeated_fields=12) -> dict:
    if not 1 <= repeats <= 20 or not 1 <= repeated_fields <= 1000:
        raise InputError("Benchmark: repeats 1..20 y repeated_fields 1..1000.")
    fields = [replace(events[0], uid=f"benchmark-{i}") for i in range(repeated_fields)]
    surface = sample_patient(equipment.patient)
    cases = []
    maps = []
    for cached in (False, True):
        durations = []
        # Warm both kernels and the spatial index; same surface/fields/reporting.
        calculate(fields, equipment, stage, surface=surface, detailed=False, cache_geometry=cached)
        for _ in range(repeats):
            start = time.perf_counter()
            calc = calculate(fields, equipment, stage, surface=surface, detailed=False, cache_geometry=cached)
            durations.append(time.perf_counter()-start)
        tracemalloc.start()
        calculate(fields, equipment, stage, surface=surface, detailed=False, cache_geometry=cached)
        _, peak = tracemalloc.get_traced_memory(); tracemalloc.stop()
        maps.append(calc.total_gy.copy())
        cases.append({"geometry_cache": cached, "times_seconds": durations,
                      "median_seconds": float(np.median(durations)),
                      "python_allocation_peak_bytes_separate_run": peak,
                      "performance": calc.summary["performance"],
                      "per_field": [ev["performance"] for ev in calc.summary["events"]]})
    return {"workload": "one beam geometry repeated; fixed prepared surface; minimal summaries",
            "repeated_fields": repeated_fields, "repeats": repeats, "cases": cases,
            "speedup_median": cases[0]["median_seconds"]/cases[1]["median_seconds"],
            "maximum_map_difference_gy": float(np.max(np.abs(maps[0]-maps[1]))),
            "memory_note": "tracemalloc mide asignaciones rastreadas en ejecución separada; no es RSS ni toda la memoria nativa.",
            "runtime": runtime_provenance()}
