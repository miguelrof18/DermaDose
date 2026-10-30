from dataclasses import replace
from copy import deepcopy
from importlib.resources import files
import numpy as np
import pytest
from skindose_verifier.config import Equipment
from skindose_verifier.synthetic import make_synthetic_rdsr
from skindose_verifier.rdsr import read_rdsr
from skindose_verifier.engine import calculate
from skindose_verifier.geometry import event_slices,rot_x,rot_y

@pytest.fixture
def model(tmp_path):
 e=Equipment.load(files('skindose_verifier').joinpath('data/equipo_sintetico.json'))
 events=read_rdsr(make_synthetic_rdsr(tmp_path/'case.dcm')).events
 return e,events

@pytest.mark.parametrize('scale',[.1,2.,10.])
def test_kerma_homogeneity(model,scale):
 eq,ev=model
 a=calculate(ev,eq,'skin_dose',detailed=False).total_gy
 b=calculate([x.scaled(scale) for x in ev],eq,'skin_dose',detailed=False).total_gy
 np.testing.assert_allclose(b,scale*a,rtol=1e-12,atol=1e-15)

def test_superposition(model):
 eq,ev=model
 total=calculate(ev,eq,'skin_dose',detailed=False).total_gy
 independent=sum(calculate([x],eq,'skin_dose',detailed=False).total_gy for x in ev)
 np.testing.assert_allclose(total,independent,rtol=1e-12,atol=1e-15)

@pytest.mark.parametrize('n',[2,7,31])
def test_sweep_conserves_reference_kerma(model,n):
 eq,ev=model;d=deepcopy(eq.raw);d.update(allow_uniform_sweep=True,sweep_samples=n)
 x=replace(ev[0],primary_end_deg=ev[0].primary_deg+40)
 parts=event_slices(x,Equipment.from_dict(d))
 assert sum(p.ka_r_gy for p in parts)==pytest.approx(x.ka_r_gy,rel=1e-14)
 assert len(parts)==n

def test_rotations_preserve_lengths_and_orientation():
 for a in [-90,0,33,180]:
  r=rot_y(a)@rot_x(37)
  np.testing.assert_allclose(r.T@r,np.eye(3),atol=1e-14)
  assert np.linalg.det(r)==pytest.approx(1.)

def test_joint_translation_preserves_free_air_map(model):
 eq,ev=model;d=deepcopy(eq.raw);shift=np.array([23.,-17.,41.])
 d['geometry']['isocenter_mm']=(np.asarray(d['geometry'].get('isocenter_mm',[0,0,0]))+shift).tolist()
 d['patient']['origin_on_table_mm']=(np.asarray(d['patient'].get('origin_on_table_mm',[0,0,0]))+shift).tolist()
 a=calculate(ev,eq,'free_air',detailed=False).total_gy
 b=calculate(ev,Equipment.from_dict(d),'free_air',detailed=False).total_gy
 np.testing.assert_allclose(a,b,rtol=1e-12,atol=1e-15)
