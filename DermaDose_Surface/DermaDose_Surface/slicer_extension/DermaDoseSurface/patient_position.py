"""Patient pose evidence for Slicer; never equate a CT origin with an isocentre."""
from hashlib import sha256
import numpy as np

POSITIONS = {
    'HFS': 'boca arriba, cabeza primero', 'FFS': 'boca arriba, pies primero',
    'HFP': 'boca abajo, cabeza primero', 'FFP': 'boca abajo, pies primero',
    'HFDR': 'sobre el lado derecho, cabeza primero',
    'FFDR': 'sobre el lado derecho, pies primero',
    'HFDL': 'sobre el lado izquierdo, cabeza primero',
    'FFDL': 'sobre el lado izquierdo, pies primero',
}
POSTURES = {'40199007': 'S', '1240000': 'P', '102535000': 'DR', '102536004': 'DL'}
LEGACY_POSTURES = {'F-10340': 'S', 'F-10310': 'P', 'F-10317': 'DR', 'F-10319': 'DL'}
ENTRANCES = {'102540008': 'HF', '102541007': 'FF'}
LEGACY_ENTRANCES = {'F-10470': 'HF', 'F-10480': 'FF'}


def concept(node):
    seq = getattr(node, 'ConceptNameCodeSequence', [])
    return (str(seq[0].CodingSchemeDesignator), str(seq[0].CodeValue)) if seq else None


def walk(node, skip_events=False):
    for child in getattr(node, 'ContentSequence', []):
        if skip_events and concept(child) == ('DCM', '113706'):
            continue
        yield child
        yield from walk(child, skip_events)


def coded_values(event, name):
    result = []
    for node in walk(event, True):
        if concept(node) == ('DCM', name):
            seq = getattr(node, 'ConceptCodeSequence', [])
            if getattr(node, 'ValueType', '') == 'CODE' and len(seq) == 1:
                result.append(dict(scheme=str(seq[0].CodingSchemeDesignator),
                                   code=str(seq[0].CodeValue)))
            else:
                result.append(dict(scheme='', code='', invalid=True))
    return result


def unique_value(values, modern, legacy):
    decoded = []
    for value in values:
        mapping = modern if value['scheme'] == 'SCT' else (
            legacy if value['scheme'] in ('SRT', 'SNM3') else {})
        decoded.append(mapping.get(value['code']))
    return decoded[0] if decoded and None not in decoded and len(set(decoded)) == 1 else None


def event_pose(event):
    evidence = dict(orientation=coded_values(event, '113743'),
                    modifier=coded_values(event, '113744'),
                    patient_table_relationship=coded_values(event, '113745'))
    posture = unique_value(evidence['modifier'], POSTURES, LEGACY_POSTURES)
    entrance = unique_value(evidence['patient_table_relationship'], ENTRANCES, LEGACY_ENTRANCES)
    recumbent = {('SCT','102538003'),('SRT','F-10450'),('SNM3','F-10450')}
    orientation_ok = all((v['scheme'],v['code']) in recumbent for v in evidence['orientation'])
    position = entrance + posture if entrance and posture and orientation_ok else None
    # Only combine fully known, event-specific codes; never parse CodeMeaning.
    uid = [str(n.UID) for n in walk(event, True)
           if concept(n) == ('DCM','113769') and getattr(n,'UID',None)]
    return dict(event_uid=uid[0] if len(uid)==1 else None,
                patient_position=position, description=POSITIONS.get(position,'posición no determinada'),
                status='known_nominal_pose' if position else 'missing_or_unsupported_or_conflicting',
                evidence=evidence)


def rdsr_pose(dataset):
    if str(getattr(dataset,'SOPClassUID','')) != '1.2.840.10008.5.1.4.1.1.88.67':
        raise ValueError('Seleccione un RDSR de rayos X.')
    events = [event_pose(node) for node in walk(dataset) if concept(node)==('DCM','113706')]
    if not events:
        raise ValueError('El RDSR no contiene eventos de irradiación.')
    positions = {e['patient_position'] for e in events}
    uniform = len(positions)==1 and None not in positions
    return dict(events=events, patient_position=next(iter(positions)) if uniform else None,
                status='uniform_nominal_pose' if uniform else 'event_specific_or_incomplete',
                root_patient_position=str(getattr(dataset,'PatientPosition','')) or None,
                registration_established=False)


def ct_geometry(ijk_to_ras, headers=()):
    matrix = np.asarray(ijk_to_ras,dtype=float)
    if matrix.shape != (4,4) or not np.isfinite(matrix).all() or not np.allclose(matrix[3],[0,0,0,1]):
        raise ValueError('Geometría del TC inválida.')
    if abs(np.linalg.det(matrix[:3,:3])) < 1e-12:
        raise ValueError('Geometría del TC singular.')
    # Slicer has already applied DICOM image orientation. Do not rotate by HFS/HFP again.
    lps = np.diag([-1.,-1.,1.,1.]) @ matrix
    positions = {str(getattr(h,'PatientPosition','')) for h in headers}
    frames = {str(getattr(h,'FrameOfReferenceUID','')) for h in headers}
    position = next(iter(positions)) if len(positions)==1 and '' not in positions else None
    def vector(header, key):
        value = getattr(header,key,None)
        return [float(v) for v in value] if value is not None else None
    return dict(ijk_to_ras_mm=matrix.tolist(), ijk_to_lps_mm=lps.tolist(),
                patient_position=position, description=POSITIONS.get(position,'posición no determinada'),
                source='DICOM' if headers else 'volume_geometry_only',
                frame_of_reference_uid=next(iter(frames)) if len(frames)==1 and '' not in frames else None,
                mixed_frames=len(frames)>1,
                dicom_image_geometry=[dict(image_orientation_patient=vector(h,'ImageOrientationPatient'),
                    image_position_patient_mm=vector(h,'ImagePositionPatient')) for h in headers],
                registration_established=False)


def read_ct_volume(volume):
    import vtk
    import slicer
    import pydicom
    matrix=vtk.vtkMatrix4x4()
    volume.GetIJKToRASMatrix(matrix)
    if volume.GetParentTransformNode():
        raise ValueError('Endurezca la transformación del TC antes de leer su geometría.')
    headers=[]
    missing=0
    for uid in (volume.GetAttribute('DICOM.instanceUIDs') or '').split():
        path=slicer.dicomDatabase.fileForInstance(uid) if slicer.dicomDatabase else ''
        if not path:
            missing+=1
            continue
        headers.append(pydicom.dcmread(path,stop_before_pixels=True,specific_tags=[
            'PatientPosition','FrameOfReferenceUID','ImageOrientationPatient','ImagePositionPatient']))
    result=ct_geometry([[matrix.GetElement(i,j) for j in range(4)] for i in range(4)],headers)
    result['unavailable_dicom_instances']=missing
    if missing:
        result['patient_position']=None
        result['description']='metadatos DICOM incompletos'
    return result


def read_rdsr_file(path):
    import pydicom
    from pathlib import Path
    result=rdsr_pose(pydicom.dcmread(path,stop_before_pixels=True))
    result['sha256']=sha256(Path(path).read_bytes()).hexdigest()
    return result


def describe(ct, rdsr=None):
    lines=['TC: '+ct['description']]
    if ct.get('mixed_frames'):
        lines.append('Aviso: hay varios marcos de referencia DICOM en los metadatos del TC.')
    if rdsr:
        position=rdsr['patient_position']
        lines.append('RDSR: '+POSITIONS.get(position,'posición variable o incompleta; revisar por evento'))
    lines.append('Registro con el equipo pendiente: hace falta localizar el isocentro en el TC y conocer los ejes del equipo.')
    lines.append('La orientación anatómica del TC ya se conserva en la malla RAS; no se gira dos veces.')
    return '\n'.join(lines)
