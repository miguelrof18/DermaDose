"""Run the reproducible v0.5 research study from the project directory.

No experimental observations are synthesized as clinical measurements.
"""
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import argparse
import json
import subprocess
import sys
import time

import numpy as np

from skindose_verifier.cli import main
from skindose_verifier.config import Equipment
from skindose_verifier.engine import calculate
from skindose_verifier.mesh_processing import mesh_quality
from skindose_verifier.patient import sample_patient
from skindose_verifier.provenance import runtime_provenance
from skindose_verifier.rdsr import read_rdsr
from skindose_verifier.uncertainty import propagate_parameters


def write(path, result):
    Path(path).write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def run_command(args):
    # Suppress the long machine-readable CLI stdout; reports remain on disk.
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()):
        status = main(args)
    if status:
        raise RuntimeError('Comando no completado: '+' '.join(args))


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--skip-tests', action='store_true')
    args = parser.parse_args()
    if not args.skip_tests:
        subprocess.run([sys.executable, '-m', 'pytest', '-q'], check=True)
    out = Path('research_results'); out.mkdir(exist_ok=True)
    print('1/6: ejemplo XCAT con aristas <=10 mm...', flush=True)
    run_command(['demo','--model','mesh','--mesh-edge-mm','10','--out','results/research_xcat','--stage','skin_dose'])
    rdsr, config = 'results/research_xcat/rdsr_sintetico.dcm', 'results/research_xcat/equipo_sintetico.json'
    eq = Equipment.load(config); events = read_rdsr(rdsr).events
    print('2/6: convergencia real 40,20,10,5 mm...', flush=True)
    run_command(['convergence',rdsr,'--config',config,'--resolutions','40,20,10,5','--out',str(out/'convergencia_xcat.json')])
    print('3/6: benchmark con y sin cache...', flush=True)
    run_command(['benchmark',rdsr,'--config',config,'--fields','12','--out',str(out/'benchmark_cache_xcat_refinado.json')])
    print('4/6: propagacion LHS con 64 muestras...', flush=True)
    run_command(['uncertainty',rdsr,'--config',config,'--spec','examples/incertidumbre_sintetica.json',
                 '--out',str(out/'uncertainty_xcat')])
    print('5/6: estabilidad con otra semilla y 256 muestras...', flush=True)
    spec = json.loads(Path('examples/incertidumbre_sintetica.json').read_text(encoding='utf-8'))
    spec.update(seed=20260929, samples=256)
    propagated = propagate_parameters(events, eq, 'skin_dose', spec)
    write(out/'uncertainty_stability_256.json', propagated.summary)
    print('6/6: escalabilidad, memoria y comparacion con malla original...', flush=True)
    records = []
    for edge in (None, 40., 20., 10.):
        data = deepcopy(eq.raw)
        data['patient'].pop('mesh_max_edge_mm', None)
        if edge is not None:
            data['patient']['mesh_max_edge_mm'] = edge
        equipment = Equipment.from_dict(data)
        begin = time.perf_counter(); surface = sample_patient(equipment.patient)
        preparation = time.perf_counter()-begin
        result = calculate(events, equipment, 'skin_dose', surface=surface)
        records.append({'max_edge_mm':edge,'preparation_seconds':preparation, 'quality':mesh_quality(surface.mesh),
                        'psd_point_gy':result.summary['psd_point_gy'], 'performance':result.summary['performance']})
    write(out/'scalability_xcat.json', {'cases':records,'runtime':runtime_provenance()})
    print('Estudio completado. Informes en research_results/.', flush=True)


if __name__ == '__main__':
    run()
