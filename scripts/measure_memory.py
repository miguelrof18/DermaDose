"""Isolated-process peak resident memory, including Python/native imports.

On Windows resource is unavailable: RSS is left null, never invented.
"""
import argparse
from copy import deepcopy
import json
import sys
import time

from skindose_verifier.config import Equipment
from skindose_verifier.engine import calculate
from skindose_verifier.rdsr import read_rdsr

p = argparse.ArgumentParser()
p.add_argument('--config', required=True)
p.add_argument('--rdsr', required=True)
p.add_argument('--edge', type=float, required=True, help='0 means original mesh')
a = p.parse_args()
cfg = deepcopy(Equipment.load(a.config).raw)
cfg['patient'].pop('mesh_max_edge_mm', None)
if a.edge:
    cfg['patient']['mesh_max_edge_mm'] = a.edge
events = read_rdsr(a.rdsr).events
started = time.perf_counter()
c = calculate(events, Equipment.from_dict(cfg), 'skin_dose', detailed=False)
elapsed = time.perf_counter()-started
try:
    import resource
    memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    memory *= 1 if sys.platform == 'darwin' else 1024
except ImportError:
    memory = None
print(json.dumps({'edge_mm':a.edge or None, 'face_count':len(c.surface.points_mm),
                  'peak_process_rss_bytes':memory, 'calculation_seconds':elapsed,
                  'scope':'fresh process; imports + surface + minimal calculation; excludes export',
                  'platform':sys.platform, 'performance':c.summary['performance']}))
