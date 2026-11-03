import sys
from pathlib import Path
import numpy as np
import pytest
from pydicom.dataset import Dataset
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'DermaDose_Surface/DermaDose_Surface/slicer_extension/DermaDoseSurface'))
from patient_position import event_pose, rdsr_pose, ct_geometry

def node(name, value, scheme='SCT'):
    d=Dataset();d.ValueType='CODE'
    n=Dataset();n.CodingSchemeDesignator='DCM';n.CodeValue=name
    d.ConceptNameCodeSequence=[n]
    v=Dataset();v.CodingSchemeDesignator=scheme;v.CodeValue=value;v.CodeMeaning='misleading text'
    d.ConceptCodeSequence=[v]
    return d

def event(posture='40199007', entrance='102540008'):
    e=node('113706','unused');e.ValueType='CONTAINER'
    orientation=node('113743','102538003')
    orientation.ContentSequence=[node('113744',posture)]
    e.ContentSequence=[orientation,node('113745',entrance)]
    return e

@pytest.mark.parametrize('posture,expected', [('40199007','HFS'),('1240000','HFP'),('102535000','HFDR'),('102536004','HFDL')])
def test_modern_posture_nested(posture,expected):
    assert event_pose(event(posture))['patient_position']==expected

def test_legacy_and_missing_not_invented():
    e=event();e.ContentSequence=[node('113744','F-10340','SRT'),node('113745','F-10480','SRT')]
    assert event_pose(e)['patient_position']=='FFS'
    e.ContentSequence.pop()
    assert event_pose(e)['patient_position'] is None

def test_unknown_scheme_and_conflict_not_parsed_by_meaning():
    e=event();e.ContentSequence.append(node('113744','1240000'))
    assert event_pose(e)['patient_position'] is None
    e=event();e.ContentSequence[1].ConceptCodeSequence[0].CodingSchemeDesignator='PRIVATE'
    assert event_pose(e)['patient_position'] is None

def test_standing_not_interpreted_as_recumbent():
    e=event();e.ContentSequence[0].ConceptCodeSequence[0].CodeValue='102539006'
    assert event_pose(e)['patient_position'] is None

def test_events_not_merged_or_inherited():
    d=Dataset();d.SOPClassUID='1.2.840.10008.5.1.4.1.1.88.67'
    d.ContentSequence=[event(),event('1240000')];d.PatientPosition='HFS'
    result=rdsr_pose(d)
    assert result['patient_position'] is None
    assert [e['patient_position'] for e in result['events']]==['HFS','HFP']
    assert not result['registration_established']

def test_ct_affine_and_no_double_rotation():
    h=Dataset();h.PatientPosition='HFP';h.FrameOfReferenceUID='1.2.3'
    matrix=np.array([[0.,-2,0,120],[2,0,0,50],[0,0,3,-200],[0,0,0,1]])
    result=ct_geometry(matrix,[h])
    assert np.allclose(result['ijk_to_ras_mm'],matrix)
    assert np.allclose(result['ijk_to_lps_mm'],np.diag([-1,-1,1,1])@matrix)
    assert not result['registration_established']
    assert ct_geometry(matrix)['patient_position'] is None
    with pytest.raises(ValueError):
        ct_geometry(np.zeros((4,4)))
