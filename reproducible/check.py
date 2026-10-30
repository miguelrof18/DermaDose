from pathlib import Path
import hashlib,json,tempfile
import numpy as np
from skindose_verifier.cli import run
ROOT=Path(__file__).resolve().parent
m=json.loads((ROOT/'manifest.json').read_text())
from importlib.resources import files
base=files('skindose_verifier').joinpath('data/benmakhlouf_2011')
for name,digest in m['factor_files_sha256'].items():
 assert hashlib.sha256(base.joinpath(name).read_bytes()).hexdigest()==digest,name
for name,digest in m['sha256'].items():
 assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
with tempfile.TemporaryDirectory() as temp:
 r=run(ROOT/'input.dcm',ROOT/'equipment.json','skin_dose',Path(temp))
 assert np.isclose(r['psd_point_gy'],m['expected_psd_gy'],rtol=1e-12,atol=1e-15)
 with np.load(r['files']['numeric_map']) as actual,np.load(ROOT/'expected_map.npz') as expected:
  np.testing.assert_allclose(actual['dose_gy'],expected['dose_gy'],rtol=1e-12,atol=1e-15)
 print('Caso reproducible correcto:',r['psd_point_gy'],'Gy')
