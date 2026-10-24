from copy import deepcopy
from types import SimpleNamespace
from pathlib import Path
import numpy as np
import pytest
from pydicom.dataset import Dataset
from skindose_verifier.rdsr import filter_stack
from skindose_verifier.domain import InputError
from skindose_verifier.surface_factors import EventSurfaceFactors,benmakhlouf_2011
from skindose_verifier.azurion import prepare_events


def node(code, value=None, scheme='DCM'):
    d=Dataset();c=Dataset();c.CodeValue=code;c.CodingSchemeDesignator='DCM';d.ConceptNameCodeSequence=[c]
    if value is not None:
        c=Dataset();c.CodeValue=value;c.CodingSchemeDesignator=scheme;d.ConceptCodeSequence=[c]
    return d


def no_filter():
    e=Dataset();c=node('113771');c.ContentSequence=[node('113772','111609')];e.ContentSequence=[c];return e


def test_explicit_no_filter_and_missing_are_distinct():
    assert filter_stack(no_filter(),'test')==(0.,[],'explicit_no_filter')
    assert filter_stack(Dataset(),'test')[0] is None


def test_contradictory_filter_is_rejected():
    e=no_filter();e.ContentSequence[0].ContentSequence.append(node('113757','C-127F9','SRT'))
    with pytest.raises(InputError,match='contradice'):filter_stack(e,'test')


def spec():
    return {'provider':'benmakhlouf_2011','inherent_al_mm':3.1,'target':'W',
            'aluminium_policy':'fixed_plus_event','field_shape_policy':'square_only',
            'out_of_domain':'error','incident_kerma_basis':'free_in_air_without_backscatter',
            'quality_basis':'event_kvp_filtration_proxy_at_surface','provenance':'test only',
            'allow_water_phantom_approximation':True}


def event(cu=.1,al=1):
    return SimpleNamespace(uid='test',kvp=80.,filtration_mm_cu=cu,field_width_mm=100.,field_height_mm=100.,
          quality={'filtration_mm_cu':'copper_with_aluminium'}, original={'filter_stack':[
          {'material_scheme':'SRT','material_code':'C-120F9','selected_mm':al}]})


def test_select_total_al_matches_published_family():
    d=benmakhlouf_2011();p=EventSurfaceFactors(spec(),d)
    actual=p.evaluate_event(event(),np.array([100.]))
    expected=d.evaluate('clinical',inherent_al_mm=4.1,kvp=80.,added_cu_mm=.1,field_side_cm=np.array([10.]))
    np.testing.assert_allclose(actual.k_bs,expected.k_bs)
    assert actual.audit['aluminium_mapping']['selected_table_al_mm']==4.1


def test_no_extrapolation_or_legacy_double_count():
    p=EventSurfaceFactors(spec(),benmakhlouf_2011())
    with pytest.raises(InputError,match='fuera del dominio'):p.evaluate_event(event(.4),np.array([100.]))
    s=spec();s.pop('aluminium_policy');p=EventSurfaceFactors(s,benmakhlouf_2011())
    with pytest.raises(InputError,match='adicional'):p.evaluate_event(event(),np.array([100.]))


def test_explicit_zero_cu_selects_unfiltered_family():
    p=EventSurfaceFactors(spec(),benmakhlouf_2011());e=event(0.,0.)
    e.quality={'filtration_mm_cu':'explicit_no_filter'}
    assert p.evaluate_event(e,np.array([100.])).audit['aluminium_mapping']['selected_table_al_mm']==3.1


def test_profile_requires_characterization_by_default():
    eq=SimpleNamespace(raw={'azurion_profile':{'name':'philips_azurion'}})
    with pytest.raises(InputError,match='caracterización'):prepare_events([],eq)
