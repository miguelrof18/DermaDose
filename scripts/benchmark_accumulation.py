"""Compare the removed stack allocation with ordered in-place accumulation.

Run from the repository: python scripts/benchmark_accumulation.py
Numbers describe aggregation only, not CT extraction or total clinical runtime.
"""
import json
import platform
import time
import tracemalloc
import numpy as np


def main(events=200, points=150_000, repeats=5):
    maps = [np.full(points, .001+i*.000001) for i in range(events)]
    def old():
        return np.sum(np.stack(maps), axis=0)
    def new():
        total = np.zeros_like(maps[0])
        for event in maps:
            total += event
        return total
    expected = old()
    cases = []
    for name, kernel in [('stack', old), ('ordered_in_place', new)]:
        kernel()  # warm-up outside timings and allocation measurement
        durations = []
        for _ in range(repeats):
            start = time.perf_counter()
            result = kernel()
            durations.append(time.perf_counter()-start)
        tracemalloc.start()
        kernel()
        peak = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
        np.testing.assert_array_equal(result, expected)
        cases.append(dict(method=name, seconds=durations, median_seconds=float(np.median(durations)),
                          additional_traced_peak_bytes=peak, max_difference_gy=float(abs(result-expected).max())))
    print(json.dumps(dict(events=events, points=points, python=platform.python_version(), numpy=np.__version__,
                          cases=cases, scope='aggregation only; retained event maps excluded from traced allocation'), indent=2))

if __name__ == '__main__':
    main()
