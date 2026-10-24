import numpy as np
import pytest
from skindose_verifier.spectral import integrate_spectrum,generate_spectrum,build_dataset,benchmark_reference
from skindose_verifier.surface_factors import benmakhlouf_2011,WaterSurfaceDataset,EventSurfaceFactors
from skindose_verifier.domain import InputError
from types import SimpleNamespace


def test_mono_node_and_normalization():
    d=benmakhlouf_2011();ref=d.evaluate('monoenergetic',energy_kev=50.,field_side_cm=10.)
    result=integrate_spectrum([50.],[20.],[.03],10.)
    assert result['combined_factor']==pytest.approx(float(ref.k_bs*ref.k_f))
    other=integrate_spectrum([50.],[20000.],[.03],10.)
    assert result['combined_factor']==other['combined_factor']


def test_product_weighting_is_not_product_of_means():
    d=benmakhlouf_2011();energies=np.array([30.,100.]);mu=np.array([.1,.02]);phi=np.array([2.,3.])
    w=energies*mu*phi;w/=w.sum()
    r=d.evaluate('monoenergetic',energy_kev=energies,field_side_cm=10.)
    result=integrate_spectrum(energies,phi,mu,10.)
    assert result['combined_factor']==pytest.approx(float(np.dot(w,r.k_bs*r.k_f)))
    assert result['k_bs']*result['k_f']==pytest.approx(result['combined_factor'])


def test_energy_midpoint_interpolates_product_on_nodes():
    d=benmakhlouf_2011();r=d.evaluate('monoenergetic',energy_kev=np.array([30.,35.]),field_side_cm=10.)
    result=integrate_spectrum([32.5],[1.],[1.],10.)
    assert result['combined_factor']==pytest.approx(float(np.mean(r.k_bs*r.k_f)))


def test_outside_domain_and_invalid_inputs():
    with pytest.raises(InputError):integrate_spectrum([4.],[1.],[1.],10.)
    with pytest.raises(InputError):integrate_spectrum([50.],[1.],[1.],4.)
    with pytest.raises(InputError):integrate_spectrum([50.],[-1.],[1.],10.)
    with pytest.raises(InputError):integrate_spectrum([50.],[float('nan')],[1.],10.)
    small=integrate_spectrum([4.,50.],[1e-10,1.],[1.,1.],10.)
    assert 0<small['excluded_kerma_fraction']<1e-4


def test_reference_agreement_real_spekpy():
    pytest.importorskip('spekpy')
    report=benchmark_reference()
    assert all(abs(x['combined_difference_percent'])<1 for x in report['cases'])


def test_generated_dataset_real_provider_and_hvl(tmp_path):
    pytest.importorskip('spekpy')
    path=tmp_path/'derived'
    meta=build_dataset(path,kvps=[80.],al_families=[4.1],copper=[.4],fields=[10.])
    assert not meta['hvl_comparison_complete']
    dataset=WaterSurfaceDataset(path)
    config={'provider':'tabulated_water','inherent_al_mm':4.1,'target':'W','field_shape_policy':'square_only',
       'out_of_domain':'error','incident_kerma_basis':'free_in_air_without_backscatter',
       'quality_basis':'event_kvp_filtration_proxy_at_surface','provenance':'test','allow_water_phantom_approximation':True}
    event=SimpleNamespace(uid='test',kvp=80.,filtration_mm_cu=.4,field_width_mm=100.,field_height_mm=100.,quality={},original={})
    pair=EventSurfaceFactors(config,dataset).evaluate_event(event,np.array([100.]))
    assert float(pair.k_bs*pair.k_f)>1
    hvl=meta['spectra'][0]['hvl1_mm_al']
    measured=[{'kvp':80.,'total_al_mm':4.1,'cu_mm':.4,'hvl1_mm_al':hvl}]
    checked=build_dataset(tmp_path/'checked',kvps=[80.],al_families=[4.1],copper=[.4],fields=[10.],measured_hvl=measured)
    assert checked['hvl_comparison_complete']
    measured[0]['hvl1_mm_al']=hvl*2
    with pytest.raises(InputError,match='HVL'):build_dataset(tmp_path/'bad',kvps=[80.],al_families=[4.1],copper=[.4],fields=[10.],measured_hvl=measured)
    assert not (tmp_path/'bad').exists()
    with pytest.raises(InputError):build_dataset(path,kvps=[80.],al_families=[4.1],copper=[.4],fields=[10.])
