"""Independent Slicer CT workflow. Run with runpy.run_path(..., run_name='__main__')."""
import json
import sys
import tempfile
import traceback
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
VERSION = '2.0.0-standalone'


def segment_body(values, threshold=-350., fill_slices=True, progress=print):
    from scipy import ndimage
    if values.ndim != 3 or not np.isfinite(threshold):
        raise ValueError('Se requiere un TC escalar 3D y un umbral finito.')
    if values.size > 160_000_000:
        raise ValueError('TC demasiado grande para este perfil. Recorte una región con Crop Volume antes de segmentar.')
    progress('Umbral de tejido y búsqueda del cuerpo principal')
    labels, count = ndimage.label(np.isfinite(values) & (values >= threshold))
    sizes = np.zeros(count+1, dtype=np.int64)
    flat = labels.reshape(-1)
    for start in range(0, flat.size, 1000000):
        sizes += np.bincount(flat[start:start+1000000], minlength=count+1)
    sizes[0] = 0
    if not sizes.any():
        raise ValueError('El umbral no selecciona tejido. Compruebe que el volumen está en HU.')
    selected = int(sizes.argmax())
    removed = int(sizes.sum()-sizes[selected])
    mask = labels == selected
    del flat, labels
    before = int(mask.sum())
    progress('Rellenando cavidades internas del contorno corporal')
    mask = ndimage.binary_fill_holes(mask)
    if fill_slices:
        # Explicit external-contour operation, not a mesh repair or morphological closing.
        # K slices are the volume's acquisition grid, not necessarily anatomical axial.
        for k in range(mask.shape[0]):
            mask[k] = ndimage.binary_fill_holes(mask[k])
    touches = [name for name,face in zip(('K−','K+','J−','J+','I−','I+'),
        (mask[0],mask[-1],mask[:,0],mask[:,-1],mask[:,:,0],mask[:,:,-1])) if face.any()]
    audit = dict(threshold_hu=float(threshold), input_shape_kji=list(values.shape),
                 connectivity=6, discarded_components=max(0,count-1),
                 discarded_voxels=removed, selected_voxels_before_fill=before,
                 filled_voxels=int(mask.sum())-before, fill_in_acquisition_slices=bool(fill_slices),
                 touches_image_faces=touches,
                 warnings=['Revise la mesa y objetos en contacto con el cuerpo.',
                    'El componente mayor puede excluir partes separadas, por ejemplo brazos.',
                    'Los límites del TC se cierran artificialmente al construir la superficie; no son piel real.'] )
    return mask.astype(np.uint8), audit


def make_surface(mask, matrix):
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk
    matrix = np.asarray(matrix,dtype=float)
    if matrix.shape!=(4,4) or not np.isfinite(matrix).all() or not np.allclose(matrix[3],[0,0,0,1]) or abs(np.linalg.det(matrix[:3,:3]))<1e-12:
        raise ValueError('Matriz IJK a RAS inválida.')
    # Contour only a tight bounding box, preserving its IJK offset exactly.
    axes = [np.flatnonzero(np.any(mask, axis=other)) for other in ((1,2),(0,2),(0,1))]
    if any(len(a)==0 for a in axes):
        raise ValueError('El segmento está vacío.')
    low = np.array([a[0] for a in axes],dtype=int)
    high = np.array([a[-1]+1 for a in axes],dtype=int)
    body=np.ascontiguousarray(mask[tuple(slice(a,b) for a,b in zip(low,high))],dtype=np.uint8)
    image=vtk.vtkImageData(); image.SetDimensions(*body.shape[::-1])
    image.GetPointData().SetScalars(numpy_to_vtk(body.ravel(),deep=True))
    pad=vtk.vtkImageConstantPad(); pad.SetInputData(image)
    nx,ny,nz=body.shape[::-1];pad.SetOutputWholeExtent(-1,nx,-1,ny,-1,nz);pad.SetConstant(0)
    contour=vtk.vtkFlyingEdges3D();contour.SetInputConnection(pad.GetOutputPort())
    contour.SetValue(0,.5);contour.ComputeNormalsOff()
    translation=np.eye(4);translation[:3,3]=low[::-1]
    effective=matrix@translation
    vm=vtk.vtkMatrix4x4()
    for i in range(4):
        for j in range(4):
            vm.SetElement(i,j,float(effective[i,j]))
    transform=vtk.vtkTransform();transform.SetMatrix(vm)
    apply=vtk.vtkTransformPolyDataFilter();apply.SetInputConnection(contour.GetOutputPort());apply.SetTransform(transform);apply.Update()
    poly=vtk.vtkPolyData();poly.DeepCopy(apply.GetOutput())
    return poly


def worker(folder):
    folder=Path(folder);job=json.loads((folder/'job.json').read_text(encoding='utf-8'))
    def progress(message): print(message,flush=True)
    values=np.load(folder/'input.npy',mmap_mode='r')
    if job['action']=='segment':
        body,audit=segment_body(values,job['threshold_hu'],job['fill_slices'],progress)
        np.save(folder/'mask.npy',body)
        (folder/'result.json').write_text(json.dumps(audit),encoding='utf-8')
    else:
        from surface_processing import simplify_surface, write_stl_checked
        import vtk
        progress('Construyendo superficie: región corporal recortada, RAS/mm')
        original=make_surface(values,job['ijk_to_ras_mm'])
        progress('Preparando malla próxima a 150 000 caras')
        poly,audit=simplify_surface(original,150000,progress,max_deviation_mm=1.,optimization_mode='fast')
        audit.update(module_version=VERSION, coordinate_system='RAS',length_unit='mm',
                     ijk_to_ras_mm=job['ijk_to_ras_mm'],registration_established=False)
        # Preserve a preview even when simplification cannot meet its constraints.
        writer=vtk.vtkXMLPolyDataWriter();writer.SetFileName(str(folder/'preview.vtp'))
        writer.SetInputData(poly);writer.SetCompressorTypeToNone();writer.EncodeAppendedDataOff()
        if not writer.Write() or writer.GetErrorCode(): raise OSError('Fallo al preparar vista 3D.')
        if audit['export_allowed']:
            audit['stl_quality']=write_stl_checked(poly,folder/'surface.stl')
        (folder/'result.json').write_text(json.dumps(audit),encoding='utf-8')
    progress('Proceso terminado')


def launch():
    import qt
    import slicer
    import vtk
    import hashlib
    import shutil
    class Window(qt.QDialog):
        def __init__(self):
            super().__init__(slicer.util.mainWindow())
            self.setWindowTitle('DermaDose CT independiente '+VERSION)
            self.setMinimumWidth(550)
            form=qt.QFormLayout(self)
            self.volume=slicer.qMRMLNodeComboBox();self.volume.nodeTypes=['vtkMRMLScalarVolumeNode']
            self.volume.noneEnabled=True;self.volume.addEnabled=False;self.volume.removeEnabled=False
            self.volume.setMRMLScene(slicer.mrmlScene);form.addRow('TC en HU:',self.volume)
            self.threshold=qt.QDoubleSpinBox();self.threshold.setRange(-1500,1000);self.threshold.value=-350
            form.addRow('Umbral inicial (HU):',self.threshold)
            self.fill=qt.QCheckBox('Rellenar aire interno también por corte del volumen');self.fill.checked=True;form.addRow(self.fill)
            self.segmentButton=qt.QPushButton('1. Segmentar cuerpo');form.addRow(self.segmentButton)
            self.editButton=qt.QPushButton('2. Revisar y corregir en Segment Editor');form.addRow(self.editButton)
            self.meshButton=qt.QPushButton('3. Preparar malla 3D (~150 000 caras)');form.addRow(self.meshButton)
            self.exportButton=qt.QPushButton('4. Exportar STL + JSON + segmentación');form.addRow(self.exportButton)
            self.cancelButton=qt.QPushButton('Cancelar proceso');self.cancelButton.enabled=False;form.addRow(self.cancelButton)
            self.log=qt.QPlainTextEdit();self.log.readOnly=True;self.log.setMinimumHeight(180);form.addRow(self.log)
            self.log.plainText='Contorno corporal externo, no una capa histológica de piel. Revise mesa y límites del TC.'
            self.process=None;self.job=None;self.cache=None;self.seg=None;self.sid=None;self.source=None;self.preview=None
            self.preparedKey=None;self.segAudit={};self.cancelled=False
            self.segmentButton.connect('clicked()',lambda:self.start('segment'))
            self.meshButton.connect('clicked()',lambda:self.start('mesh'))
            self.editButton.connect('clicked()',self.edit)
            self.exportButton.connect('clicked()',self.export)
            self.cancelButton.connect('clicked()',self.cancel)
            volumes=[v for v in slicer.util.getNodesByClass('vtkMRMLScalarVolumeNode') if v.GetImageData()]
            if len(volumes)==1:self.volume.setCurrentNode(volumes[0])
        def controls(self,running):
            for widget in (self.segmentButton,self.meshButton,self.editButton,self.exportButton,self.volume,self.threshold,self.fill):widget.enabled=not running
            self.cancelButton.enabled=running
        def key(self):
            if self.seg is None or self.seg.GetScene() is None or self.source.GetScene() is None:
                raise ValueError('Segmente primero el TC.')
            if self.seg.GetParentTransformNode() or self.source.GetParentTransformNode():
                raise ValueError('Endurezca las transformaciones antes de preparar la malla.')
            segment=self.seg.GetSegmentation().GetSegment(self.sid)
            if segment is None:raise ValueError('El segmento fue eliminado.')
            binary=segment.GetRepresentation('Binary labelmap')
            return (self.seg.GetMTime(),binary.GetMTime(),self.source.GetMTime(),self.source.GetImageData().GetMTime())
        def start(self,action):
            if self.process is not None:return
            try:
                volume=self.volume.currentNode()
                if volume is None or volume.GetImageData() is None:raise ValueError('Seleccione un TC escalar en HU.')
                if volume.GetParentTransformNode():raise ValueError('Endurezca la transformación del TC primero.')
                if action=='mesh' and volume!=self.source:raise ValueError('Seleccione el TC usado para segmentar.')
                if action=='segment':values=slicer.util.arrayFromVolume(volume)
                else:values=slicer.util.arrayFromSegmentBinaryLabelmap(self.seg,self.sid,volume)
                if values.ndim!=3 or values.size>160_000_000:raise ValueError('Recorte el TC con Crop Volume: tamaño no admitido por este perfil.')
                executable=Path(slicer.app.applicationDirPath())/'PythonSlicer.exe'
                if not executable.exists():executable=executable.with_suffix('')
                if not executable.exists():raise FileNotFoundError('No se encuentra PythonSlicer junto a Slicer.')
                self.job=tempfile.TemporaryDirectory(prefix='dermadose_ct_');folder=Path(self.job.name)
                matrix=vtk.vtkMatrix4x4();volume.GetIJKToRASMatrix(matrix)
                config=dict(action=action,threshold_hu=float(self.threshold.value),fill_slices=bool(self.fill.checked),
                    ijk_to_ras_mm=[[matrix.GetElement(i,j) for j in range(4)] for i in range(4)])
                (folder/'job.json').write_text(json.dumps(config),encoding='utf-8')
                self.action=action;self.jobVolume=volume;self.volumeTime=(volume.GetMTime(),volume.GetImageData().GetMTime())
                self.jobKey=self.key() if action=='mesh' else None
                self.controls(True);self.log.plainText='Copiando datos al proceso de cálculo…';slicer.app.processEvents()
                np.save(folder/'input.npy',values)
                self.cancelled=False;self.process=qt.QProcess(self)
                self.process.setProcessChannelMode(qt.QProcess.MergedChannels)
                self.process.connect('readyReadStandardOutput()',self.read)
                self.process.connect('finished(int,QProcess::ExitStatus)',self.finished)
                self.process.connect('errorOccurred(QProcess::ProcessError)',self.error)
                self.process.start(str(executable),['-u',str(HERE/'dermadose_ct.py'),'--worker',str(folder)])
            except Exception as exc:
                self.controls(False)
                if self.job:self.job.cleanup();self.job=None
                self.log.plainText=str(exc)
        def read(self):
            if self.process:self.log.appendPlainText(bytes(self.process.readAllStandardOutput()).decode('utf-8',errors='replace'))
        def error(self,code):
            if code==qt.QProcess.FailedToStart:self.log.appendPlainText(self.process.errorString());self.finished(-1,1)
        def cancel(self):
            if self.process:self.cancelled=True;self.process.kill()
        def finished(self,code,status):
            self.read()
            keep=False
            try:
                if self.cancelled:self.log.appendPlainText('Cancelado; se conserva la segmentación anterior.');return
                if code!=0:raise RuntimeError('El proceso falló. El detalle aparece arriba; no se ha descartado la segmentación.')
                volume=self.jobVolume
                if volume.GetScene() is None or volume.GetImageData() is None or (volume.GetMTime(),volume.GetImageData().GetMTime())!=self.volumeTime:
                    raise ValueError('El TC cambió durante el proceso. Repita la operación.')
                folder=Path(self.job.name);audit=json.loads((folder/'result.json').read_text(encoding='utf-8'))
                if self.action=='segment':
                    seg=slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode','DermaDose CT body')
                    try:
                        seg.CreateDefaultDisplayNodes();seg.SetReferenceImageGeometryParameterFromVolumeNode(volume)
                        sid=seg.GetSegmentation().AddEmptySegment('Body','Contorno corporal')
                        slicer.util.updateSegmentBinaryLabelmapFromArray(np.load(folder/'mask.npy'),seg,sid,volume)
                        seg.GetDisplayNode().SetVisibility3D(False)
                    except Exception:
                        slicer.mrmlScene.RemoveNode(seg);raise
                    if self.seg and self.seg.GetScene():self.seg.GetDisplayNode().SetVisibility(False)
                    self.seg,self.sid,self.source=seg,sid,volume;self.segAudit=audit;self.preparedKey=None
                    if self.preview:self.preview.GetDisplayNode().SetVisibility(False)
                    slicer.util.setSliceViewerLayers(background=volume)
                    self.log.appendPlainText('Segmentación lista en los cortes. Revise mesa, brazos y cobertura antes de preparar 3D.\n'+'\n'.join(audit['warnings']))
                else:
                    if self.key()!=self.jobKey:raise ValueError('El segmento cambió durante el cálculo; prepare la malla de nuevo.')
                    reader=vtk.vtkXMLPolyDataReader();reader.SetFileName(str(folder/'preview.vtp'));reader.Update()
                    if not self.preview or self.preview.GetScene() is None:
                        self.preview=slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode','DermaDose CT mesh preview');self.preview.CreateDefaultDisplayNodes()
                    poly=vtk.vtkPolyData();poly.DeepCopy(reader.GetOutput());self.preview.SetAndObservePolyData(poly)
                    self.preview.GetDisplayNode().SetColor(.85,.72,.6);self.preview.GetDisplayNode().SetVisibility(True)
                    self.seg.GetDisplayNode().SetVisibility3D(False)
                    if self.cache:self.cache.cleanup()
                    self.cache=self.job;keep=True;self.preparedKey=self.key()
                    self.log.appendPlainText(f"Malla: {audit['input_triangles']:,} → {audit['output_triangles']:,} caras; desviación muestreada {audit['max_vertex_surface_distance_mm']:.3f} mm.")
                    if not audit['export_allowed']:self.log.appendPlainText('No se obtuvo una reducción válida en el margen. La segmentación sigue disponible para corregirla o guardarla desde Slicer.')
                    slicer.util.resetThreeDViews()
            except Exception as exc:self.log.appendPlainText(str(exc))
            finally:
                self.process.deleteLater();self.process=None
                if not keep:self.job.cleanup()
                self.job=None;self.controls(False)
        def edit(self):
            try:
                self.key()
                if self.preview:self.preview.GetDisplayNode().SetVisibility(False)
                slicer.util.selectModule('SegmentEditor')
                editor=slicer.modules.segmenteditor.widgetRepresentation().self().editor
                editor.setSegmentationNode(self.seg);editor.setSourceVolumeNode(self.source);editor.setCurrentSegmentID(self.sid)
                self.log.appendPlainText('Corrija con Scissors/Islands. Tras editar, vuelva a preparar la malla 3D.')
            except Exception as exc:self.log.appendPlainText(str(exc))
        def export(self, destination=None):
            try:
                if self.cache is None or self.key()!=self.preparedKey:raise ValueError('Prepare la malla después de la última edición.')
                folder=Path(self.cache.name);audit=json.loads((folder/'result.json').read_text(encoding='utf-8'))
                if not audit['export_allowed']:raise ValueError('Corrija el segmento: no hay una malla reducida válida para exportar.')
                if destination is None:
                    destination=qt.QFileDialog.getSaveFileName(self,'Exportar malla revisada','DermaDose_CT.stl','STL (*.stl)')
                if not destination:return
                target=Path(destination)
                if target.suffix.lower()!='.stl':raise ValueError('El nombre debe terminar en .stl.')
                paths=[target,target.with_suffix('.json'),target.with_suffix('.seg.nrrd')]
                if any(p.exists() for p in paths):raise FileExistsError('Ya existe alguna salida. Elija otro nombre.')
                with tempfile.TemporaryDirectory(dir=str(target.parent)) as staging:
                    staging=Path(staging);segfile=staging/'body.seg.nrrd'
                    if not slicer.util.saveNode(self.seg,str(segfile)):raise OSError('No se pudo guardar la segmentación editable.')
                    audit.update(segmentation=self.segAudit,mesh_sha256=hashlib.sha256((folder/'surface.stl').read_bytes()).hexdigest(),
                                 surface_definition='external body contour; review required',manual_edits='included from current binary labelmap')
                    (staging/'audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
                    created=[]
                    try:
                        for src,dst in zip([folder/'surface.stl',staging/'audit.json',segfile],paths):
                            with dst.open('xb') as output:
                                created.append(dst)
                                with src.open('rb') as incoming:shutil.copyfileobj(incoming,output)
                    except Exception:
                        for path in created:path.unlink(missing_ok=True)
                        raise
                self.log.appendPlainText('Guardados STL, JSON y .seg.nrrd editable en '+str(target.parent))
            except Exception as exc:self.log.appendPlainText(str(exc))
        def closeEvent(self,event):
            if self.process:
                self.process.blockSignals(True);self.process.kill();self.process.waitForFinished(3000)
            if self.job:self.job.cleanup()
            if self.cache:self.cache.cleanup()
            event.accept()
    previous=getattr(slicer,'_dermadose_ct_window',None)
    if previous:previous.close()
    slicer._dermadose_ct_window=Window();slicer._dermadose_ct_window.show()
    return slicer._dermadose_ct_window


def self_test():
    """Run INSIDE real Slicer; preserves the existing scene and exercises the GUI API."""
    import time
    import slicer
    import vtk
    data=np.full((24,30,36),-1000,dtype=np.int16)
    data[4:20,5:25,6:30]=40
    data[9:14,11:18,12:22]=-1000
    volume=slicer.util.addVolumeFromArray(data,name='DermaDose self-test CT')
    volume.SetSpacing(1.5,2.,3.)
    window=launch();window.volume.setCurrentNode(volume)
    nodes=[]
    def wait():
        deadline=time.monotonic()+120
        while window.process is not None and time.monotonic()<deadline:
            slicer.app.processEvents();time.sleep(.02)
        if window.process is not None:raise TimeoutError('El proceso no termina en 120 s.')
    try:
        window.start('segment');wait()
        assert window.seg is not None, window.log.plainText
        nodes.append(window.seg)
        window.start('mesh');wait()
        assert window.cache is not None, window.log.plainText
        nodes.append(window.preview)
        with tempfile.TemporaryDirectory() as folder:
            target=Path(folder)/'test.stl'
            window.export(str(target))
            assert target.exists() and target.with_suffix('.seg.nrrd').exists(),window.log.plainText
            reader=vtk.vtkSTLReader();reader.SetFileName(str(target));reader.Update()
            assert reader.GetOutput().GetNumberOfPolys()>4
            metadata=json.loads(target.with_suffix('.json').read_text(encoding='utf-8'))
            assert metadata['coordinate_system']=='RAS'
        print('DERMADOSE_SELF_TEST_OK: segmentación, proceso externo, malla y exportación.')
    finally:
        window.close()
        for node in nodes+[volume]:
            if node is not None and node.GetScene():slicer.mrmlScene.RemoveNode(node)


if __name__=='__main__':
    if len(sys.argv)>2 and sys.argv[1]=='--worker':
        try:worker(sys.argv[2])
        except Exception:traceback.print_exc();sys.exit(1)
    else:launch()
