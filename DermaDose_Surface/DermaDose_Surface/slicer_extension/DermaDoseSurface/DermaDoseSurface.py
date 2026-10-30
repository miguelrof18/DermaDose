"""3D Slicer module: CT external contour for research use in DermaDose."""

import hashlib
import json
import tempfile
from contextlib import contextmanager
from vtk.util.numpy_support import vtk_to_numpy
from pathlib import Path

import numpy as np
import qt
import slicer
import vtk
from scipy import ndimage
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule,
    ScriptedLoadableModuleLogic,
    ScriptedLoadableModuleTest,
    ScriptedLoadableModuleWidget,
)


MODULE_VERSION = "1.3.0-150k"

def surface_quality(poly):
    connectivity = vtk.vtkPolyDataConnectivityFilter()
    connectivity.SetInputData(poly)
    connectivity.SetExtractionModeToAllRegions()
    connectivity.Update()
    counts = {}
    for name in ("boundary_edges", "nonmanifold_edges"):
        edges = vtk.vtkFeatureEdges()
        edges.SetInputData(poly)
        edges.FeatureEdgesOff()
        edges.ManifoldEdgesOff()
        edges.BoundaryEdgesOff()
        edges.NonManifoldEdgesOff()
        if name == "boundary_edges":
            edges.BoundaryEdgesOn()
        else:
            edges.NonManifoldEdgesOn()
        edges.Update()
        counts[name] = int(edges.GetOutput().GetNumberOfCells())
    triangle_filter = vtk.vtkTriangleFilter()
    triangle_filter.SetInputData(poly)
    triangle_filter.PassVertsOff()
    triangle_filter.PassLinesOff()
    triangle_filter.Update()
    triangles = triangle_filter.GetOutput()
    degenerate = 0
    finite = False
    volume = 0.0
    if triangles.GetNumberOfPoints() and triangles.GetNumberOfPolys():
        points = vtk_to_numpy(triangles.GetPoints().GetData())
        faces = vtk_to_numpy(triangles.GetPolys().GetConnectivityArray()).reshape(-1, 3)
        finite = bool(np.isfinite(points).all())
        for start in range(0, len(faces), 100000):
            vertices = points[faces[start:start + 100000]].astype(np.float64)
            cross = np.cross(vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0])
            degenerate += int(np.count_nonzero(np.linalg.norm(cross, axis=1) <= 1e-12))
            volume += float(np.einsum("ij,ij->", vertices[:, 0], cross) / 6.0)
    return {
        "finite_coordinates": finite,
        "degenerate_triangles": degenerate,
        "signed_volume_mm3": volume,
        "triangles": int(poly.GetNumberOfPolys()),
        "components": int(connectivity.GetNumberOfExtractedRegions()),
        **counts,
    }


def clean_surface(poly):
    # Only weld exactly coincident vertices; do not close holes or move anatomy.
    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(poly)
    clean.PointMergingOn()
    clean.ToleranceIsAbsoluteOn()
    clean.SetAbsoluteTolerance(0.0)
    clean.ConvertPolysToLinesOff()
    clean.ConvertLinesToPointsOff()
    clean.Update()
    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(clean.GetOutputPort())
    triangles.PassVertsOff()
    triangles.PassLinesOff()
    triangles.Update()
    output = vtk.vtkPolyData()
    output.DeepCopy(triangles.GetOutput())
    return output


def valid_surface(quality):
    return (quality["triangles"] >= 4 and quality["components"] == 1
            and quality["boundary_edges"] == 0 and quality["nonmanifold_edges"] == 0
            and quality["finite_coordinates"] and quality["degenerate_triangles"] == 0
            and abs(quality["signed_volume_mm3"]) > 1e-9)


def simplify_surface(poly, target_triangles=150000):
    """Try conservative reductions; retain the original if all reductions fail."""
    if target_triangles < 100:
        raise ValueError("El objetivo mínimo es de 100 triángulos.")
    source = clean_surface(poly)
    before = surface_quality(source)
    if not valid_surface(before):
        raise ValueError("Superficie original inválida: " + str(before)
                         + ". Corrija la segmentación antes de simplificar.")
    count = before["triangles"]
    result = source
    attempts = []
    reduced = False
    def orient(candidate):
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputData(candidate)
        normals.SplittingOff()
        normals.ConsistencyOn()
        normals.AutoOrientNormalsOn()
        normals.Update()
        return clean_surface(normals.GetOutput())

    budgets = []
    if count > target_triangles:
        budgets = sorted(set([int(target_triangles), int((target_triangles * count) ** 0.5),
                              max(int(target_triangles), int(count * 0.75))]))
    chosen_method = "unchanged"
    found = False
    for budget in budgets:
        if budget >= count:
            continue
        # Try both available VTK decimators; accept only independently checked outputs.
        for method in ("vtkDecimatePro", "vtkQuadricDecimation"):
            decimator = getattr(vtk, method)()
            decimator.SetInputData(source)
            decimator.SetTargetReduction(1.0 - budget / count)
            if method == "vtkDecimatePro":
                decimator.PreserveTopologyOn()
                decimator.SplittingOff()
                decimator.BoundaryVertexDeletionOff()
            decimator.Update()
            candidate = orient(clean_surface(decimator.GetOutput()))
            quality = surface_quality(candidate)
            attempts.append({"method": method, "requested_triangles": budget, **quality})
            if valid_surface(quality) and quality["triangles"] < count:
                result, chosen_method, found = candidate, method, True
                break
        if found:
            break
    if not found:
        result = orient(source)
    after = surface_quality(result)
    if not valid_surface(after):
        raise ValueError("Superficie original: fallo al orientar normales. Diagnóstico: " + str(after))
    return result, {
        "method": chosen_method,
        "target_triangles": int(target_triangles),
        "input_triangles": count, "output_triangles": after["triangles"],
        "target_reached": after["triangles"] <= target_triangles,
        "fallback_original": count > target_triangles and not found,
        "attempts": attempts,
    }


def write_stl_checked(poly, path):
    """Validate the bytes actually written, before exposing an STL to the user."""
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(str(path))
    writer.SetFileTypeToBinary()
    writer.SetHeader("DermaDose SPACE=RAS; coordinates in mm")
    writer.SetInputData(poly)
    if not writer.Write() or writer.GetErrorCode() != 0:
        raise OSError("No se pudo escribir el STL.")
    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(path))
    reader.MergingOn()
    reader.Update()
    quality = surface_quality(reader.GetOutput())
    if not valid_surface(quality):
        raise ValueError("El STL guardado pierde integridad al convertir a precisión STL: " + str(quality))
    return quality


class DermaDoseSurface(ScriptedLoadableModule):
    def __init__(self, parent):
        super().__init__(parent)
        parent.title = "DermaDose Surface"
        parent.categories = ["DermaDose"]
        parent.contributors = ["DermaDose"]
        parent.helpText = (
            "Extract an external CT contour, inspect it, and export a RAS/mm STL "
            f"with an audit JSON. Version {MODULE_VERSION}. This is not a skin-layer segmentation."
        )


class DermaDoseSurfaceWidget(ScriptedLoadableModuleWidget):
    def setup(self):
        super().setup()
        self.logic = DermaDoseSurfaceLogic()
        form = qt.QFormLayout()
        self.layout.addLayout(form)

        self.volume = slicer.qMRMLNodeComboBox()
        self.volume.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.volume.noneEnabled = True
        self.volume.addEnabled = False
        self.volume.removeEnabled = False
        self.volume.setMRMLScene(slicer.mrmlScene)
        form.addRow("TC (volumen):", self.volume)
        form.addRow(qt.QLabel(f"DermaDose Surface {MODULE_VERSION}"))
        modulePath = qt.QLineEdit(str(Path(__file__).resolve()))
        modulePath.readOnly = True
        form.addRow("Archivo cargado:", modulePath)

        self.threshold = qt.QDoubleSpinBox()
        self.threshold.minimum = -1500
        self.threshold.maximum = 3000
        self.threshold.value = -350
        self.threshold.suffix = " HU"
        self.threshold.toolTip = "Umbral inicial orientativo; revise el contorno en los cortes."
        form.addRow("Umbral mínimo:", self.threshold)

        self.targetTriangles = qt.QSpinBox()
        self.targetTriangles.minimum = 100
        self.targetTriangles.maximum = 10000000
        self.targetTriangles.singleStep = 10000
        self.targetTriangles.value = 150000
        self.targetTriangles.toolTip = "Objetivo aproximado; se conserva la topología aunque se supere el objetivo."
        self.targetTriangles.connect("valueChanged(int)", self.onTargetChanged)
        form.addRow("Objetivo de triángulos:", self.targetTriangles)

        self.extractButton = qt.QPushButton("Extraer superficie")
        self.extractButton.connect("clicked()", self.onExtract)
        form.addRow(self.extractButton)

        self.editButton = qt.QPushButton("Abrir Segment Editor para corregir")
        self.editButton.connect("clicked()", self.onEdit)
        form.addRow(self.editButton)

        self.qaButton = qt.QPushButton("Revisar malla")
        self.qaButton.connect("clicked()", self.onCheck)
        form.addRow(self.qaButton)

        self.report = qt.QPlainTextEdit()
        self.report.readOnly = True
        self.report.setMinimumHeight(130)
        form.addRow("Comprobación:", self.report)

        self.acknowledge = qt.QCheckBox(
            "He revisado contorno, mesa y cobertura del TC; acepto los avisos."
        )
        form.addRow(self.acknowledge)

        self.exportButton = qt.QPushButton("Exportar STL y metadatos")
        self.exportButton.connect("clicked()", self.onExport)
        form.addRow(self.exportButton)
        self.layout.addStretch(1)
        self._reviewedKey = None

    @contextmanager
    def busy(self):
        controls = (self.extractButton, self.editButton, self.qaButton, self.exportButton,
                    self.targetTriangles, self.volume, self.threshold, self.acknowledge)
        self._reviewedKey = None
        self.acknowledge.checked = False
        for control in controls:
            control.setEnabled(False)
        slicer.util.mainWindow().setCursor(qt.Qt.WaitCursor)
        slicer.app.processEvents()
        try:
            yield
        finally:
            slicer.util.mainWindow().unsetCursor()
            for control in controls:
                control.setEnabled(True)

    def onTargetChanged(self, value):
        self.logic.targetTriangles = int(value)
        self.acknowledge.checked = False
        self.report.plainText = "Objetivo cambiado. Pulse Revisar malla para actualizar la vista previa."

    def onExtract(self):
        try:
            with self.busy():
                self.logic.extract(self.volume.currentNode(), self.threshold.value)
            self.acknowledge.checked = False
            self.onCheck()
            slicer.util.selectModule("DermaDoseSurface")
        except Exception as exc:
            slicer.util.errorDisplay(str(exc))

    def onEdit(self):
        if self.logic.segmentationNode is None:
            slicer.util.errorDisplay("Extraiga primero el contorno.")
            return
        if self.logic.previewNode is not None:
            self.logic.previewNode.GetDisplayNode().SetVisibility(False)
        self.logic.segmentationNode.GetDisplayNode().SetVisibility3D(True)
        slicer.util.selectModule("SegmentEditor")
        editor = slicer.modules.segmenteditor.widgetRepresentation().self().editor
        editor.setSegmentationNode(self.logic.segmentationNode)
        editor.setSourceVolumeNode(self.logic.volumeNode)
        editor.setCurrentSegmentID(self.logic.segmentId)
        self.acknowledge.checked = False

    def onCheck(self):
        try:
            with self.busy():
                result = self.logic.inspect()
                self.logic.showPreview()
            self._reviewedKey = self.logic._preparedKey
            lines = [
                f"Triángulos originales: {result['original_triangles']}",
                f"Triángulos de salida: {result['triangles']} (objetivo {self.logic.targetTriangles})",
                f"Componentes de superficie: {result['components']}",
                f"Aristas no manifold: {result['nonmanifold_edges']}",
                f"Aristas abiertas: {result['boundary_edges']}",
                f"Triángulos degenerados: {result['degenerate_triangles']}",
                f"Método: {self.logic.decimation['method']}",
                *["AVISO: " + warning for warning in result["warnings"]],
            ]
            self.report.plainText = "\n".join(lines)
            self.acknowledge.checked = False
        except Exception as exc:
            slicer.util.errorDisplay(str(exc))

    def onExport(self):
        if not self.acknowledge.checked:
            slicer.util.errorDisplay("Revise la malla y marque la casilla antes de exportar.")
            return
        target = qt.QFileDialog.getSaveFileName(
            slicer.util.mainWindow(), "Exportar superficie", "DermaDose_surface.stl", "STL (*.stl)"
        )
        if not target:
            return
        try:
            self.logic.preparedSurface()
            if self._reviewedKey != self.logic._preparedKey:
                self.acknowledge.checked = False
                raise ValueError("La malla cambió desde la revisión. Pulse Revisar malla y compruebe la nueva vista previa.")
            with self.busy():
                self.logic.export(Path(target))
            slicer.util.infoDisplay(f"Exportados {target} y su archivo .json")
        except Exception as exc:
            slicer.util.errorDisplay(str(exc))


class DermaDoseSurfaceLogic(ScriptedLoadableModuleLogic):
    def __init__(self):
        super().__init__()
        self.segmentationNode = None
        self.volumeNode = None
        self.segmentId = None
        self.thresholdHU = None
        self.targetTriangles = 150000
        self.previewNode = None
        self._preparedKey = None
        self._preparedPoly = None
        self.decimation = None

    def extract(self, volume, threshold_hu):
        if volume is None or volume.GetImageData() is None:
            raise ValueError("Cargue y seleccione una serie de TC.")
        if volume.GetParentTransformNode():
            raise ValueError("Endurezca la transformación del TC antes de extraer la malla.")
        values = slicer.util.arrayFromVolume(volume)
        if values.ndim != 3:
            raise ValueError("El volumen seleccionado debe ser escalar y tridimensional.")
        mask = values >= threshold_hu
        if not np.any(mask):
            raise ValueError("El umbral no selecciona ningún vóxel.")
        # Keep only the largest 3D connected body. Table contact still needs manual review.
        labels, _ = ndimage.label(mask, structure=ndimage.generate_binary_structure(3, 1))
        sizes = np.bincount(labels.ravel())
        sizes[0] = 0
        body = labels == int(np.argmax(sizes))
        del labels, mask
        body = ndimage.binary_fill_holes(body).astype(np.uint8)
        previous = self.segmentationNode
        seg = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", "DermaDose external contour")
        try:
            seg.CreateDefaultDisplayNodes()
            seg.SetReferenceImageGeometryParameterFromVolumeNode(volume)
            segment_id = seg.GetSegmentation().AddEmptySegment("ExternalContour", "External contour")
            slicer.util.updateSegmentBinaryLabelmapFromArray(body, seg, segment_id, volume)
            if not seg.CreateClosedSurfaceRepresentation():
                raise RuntimeError("No se pudo generar la superficie cerrada.")
            seg.GetDisplayNode().SetAllSegmentsVisibility(True)
            if previous is not None and previous.GetScene() is not None:
                previous.GetDisplayNode().SetVisibility(False)
            self._preparedKey = None
            self._preparedPoly = None
            if self.previewNode is not None and self.previewNode.GetScene() is not None:
                self.previewNode.GetDisplayNode().SetVisibility(False)
            self.segmentationNode, self.volumeNode = seg, volume
            self.segmentId, self.thresholdHU = segment_id, float(threshold_hu)
            slicer.util.resetThreeDViews()
            return seg
        except Exception:
            slicer.mrmlScene.RemoveNode(seg)
            raise

    def _surface(self):
        if (not self.segmentationNode or not self.segmentId or not self.volumeNode
                or self.segmentationNode.GetScene() is None or self.volumeNode.GetScene() is None):
            raise ValueError("Extraiga primero la superficie del TC.")
        if self.volumeNode.GetParentTransformNode() or self.segmentationNode.GetParentTransformNode():
            raise ValueError("Endurezca las transformaciones antes de exportar.")
        self.segmentationNode.CreateClosedSurfaceRepresentation()
        segment = self.segmentationNode.GetSegmentation().GetSegment(self.segmentId)
        if segment is None:
            raise ValueError("El segmento ya no existe; vuelva a extraer la superficie.")
        poly = segment.GetRepresentation("Closed surface")
        if poly is None or poly.GetNumberOfPolys() == 0:
            raise ValueError("No se pudo obtener la superficie del segmento.")
        return poly

    def preparedSurface(self):
        # Read the live representation without copying millions of triangles on cache hits.
        representation = self._surface()
        key = (self.segmentationNode.GetID(), representation.GetMTime(),
               representation.GetPoints().GetMTime(), representation.GetPolys().GetMTime(), self.targetTriangles)
        if key != self._preparedKey:
            self._preparedPoly, self.decimation = simplify_surface(representation, self.targetTriangles)
            self._preparedKey = key
        return self._preparedPoly

    def showPreview(self):
        poly = self.preparedSurface()
        if self.previewNode is None or self.previewNode.GetScene() is None:
            self.previewNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "DermaDose export preview")
            self.previewNode.CreateDefaultDisplayNodes()
        self.previewNode.SetAndObservePolyData(poly)
        self.previewNode.GetDisplayNode().SetColor(0.85, 0.7, 0.55)
        self.previewNode.GetDisplayNode().SetVisibility(True)
        self.segmentationNode.GetDisplayNode().SetVisibility3D(False)

    def inspect(self):
        poly = self.preparedSurface()
        quality = surface_quality(poly)
        mask = slicer.util.arrayFromSegmentBinaryLabelmap(
            self.segmentationNode, self.segmentId, self.volumeNode
        )
        # Axis order is KJI. A contour touching any image face may be truncated.
        touches = [name for name, face in zip(
            ("K mínimo", "K máximo", "J mínimo", "J máximo", "I mínimo", "I máximo"),
            (mask[0], mask[-1], mask[:, 0], mask[:, -1], mask[:, :, 0], mask[:, :, -1]),
        ) if np.any(face)]
        warnings = ["La mesa puede quedar conectada al paciente; revise cortes y vista 3D."]
        if touches:
            warnings.append("El contorno alcanza el límite del TC: " + ", ".join(touches))
        if quality["components"] != 1:
            warnings.append("Hay varios componentes de superficie; revise objetos extraños.")
        if quality["boundary_edges"]:
            warnings.append("Hay aristas abiertas; revise la integridad de la malla.")
        if self.decimation.get("fallback_original"):
            warnings.append("No se obtuvo una reducción válida: la vista previa y exportación conservan la malla original. No acelerará el cálculo.")
        if not self.decimation["target_reached"]:
            warnings.append("Objetivo no alcanzado: se conserva la topología. Puede exportar con más triángulos.")
        return {
            **quality,
            "original_triangles": self.decimation["input_triangles"],
            "warnings": warnings,
        }

    def export(self, path):
        path = Path(path)
        if path.suffix.lower() != ".stl":
            raise ValueError("Seleccione un nombre terminado en .stl.")
        json_path = path.with_suffix(".json")
        if path.exists() or json_path.exists():
            raise FileExistsError("Ya existe el STL o su JSON; elija un nombre nuevo (por ejemplo superficie_150k_v3.stl).")
        qa = self.inspect()
        poly = self.preparedSurface()
        matrix = vtk.vtkMatrix4x4()
        self.volumeNode.GetIJKToRASMatrix(matrix)
        # Stage in the destination filesystem. Never overwrite an existing file.
        with tempfile.TemporaryDirectory(prefix="dermadose_", dir=str(path.parent)) as folder:
            staged = Path(folder) / "surface.stl"
            disk_quality = write_stl_checked(poly, staged)
            metadata = {
                "schema_version": 4, "module_version": MODULE_VERSION,
                "decimation": self.decimation,
                "surface_definition": "external CT contour, not histological skin",
                "mesh_file": path.name,
                "mesh_sha256": hashlib.sha256(staged.read_bytes()).hexdigest(),
                "coordinate_system": "RAS", "length_unit": "mm",
                "ijk_to_ras_mm": [[matrix.GetElement(i, j) for j in range(4)] for i in range(4)],
                "threshold_minimum_hu": self.thresholdHU,
                "processing": ["threshold", "largest_6_connected_component", "fill_enclosed_holes",
                               "checked_decimation", "STL_readback_validation"],
                "manual_edits": "Manual changes are not automatically recorded",
                "quality": qa, "stl_quality": disk_quality,
                "rdsr_registration": "Not established; verify geometry separately before dose calculation",
            }
            created_stl = created_json = False
            try:
                with path.open("xb") as output:
                    created_stl = True
                    with staged.open("rb") as source:
                        import shutil
                        shutil.copyfileobj(source, output)
                with json_path.open("x", encoding="utf-8") as output:
                    created_json = True
                    json.dump(metadata, output, indent=2, ensure_ascii=False)
                    output.write("\n")
            except Exception:
                if created_stl:
                    path.unlink(missing_ok=True)
                if created_json:
                    json_path.unlink(missing_ok=True)
                raise
        return metadata


class DermaDoseSurfaceTest(ScriptedLoadableModuleTest):
    def runTest(self):
        self.setUp()
        self.testSyntheticContour()
        self.testDecimation()

    def testSyntheticContour(self):
        data = np.full((32, 32, 32), -1000, dtype=np.int16)
        data[5:27, 6:26, 7:25] = 50
        data[11:14, 12:15, 13:16] = -1000  # enclosed cavity must be filled
        data[1:3, 1:3, 1:3] = 50  # disconnected object must be discarded
        volume = slicer.util.addVolumeFromArray(data, name="DermaDose synthetic CT")
        logic = DermaDoseSurfaceLogic()
        try:
            logic.extract(volume, -350)
            qa = logic.inspect()
            self.assertGreater(qa["triangles"], 0)
            self.assertEqual(qa["components"], 1)
            self.assertEqual(qa["boundary_edges"], 0)
            with tempfile.TemporaryDirectory() as folder:
                output = Path(folder) / "surface.stl"
                metadata = logic.export(output)
                self.assertTrue(output.exists())
                self.assertEqual(metadata["mesh_sha256"], hashlib.sha256(output.read_bytes()).hexdigest())
                self.assertEqual(metadata["coordinate_system"], "RAS")
                self.assertEqual(json.loads(output.with_suffix(".json").read_text())["length_unit"], "mm")
                with self.assertRaises(FileExistsError):
                    logic.export(output)
        finally:
            if logic.segmentationNode:
                slicer.mrmlScene.RemoveNode(logic.segmentationNode)
            slicer.mrmlScene.RemoveNode(volume)


    def testDecimation(self):
        sphere = vtk.vtkSphereSource()
        sphere.SetThetaResolution(400)
        sphere.SetPhiResolution(300)
        sphere.Update()
        reduced, info = simplify_surface(sphere.GetOutput(), 150000)
        self.assertGreater(info["input_triangles"], 150000)
        self.assertLessEqual(info["output_triangles"], 150000)
        quality = surface_quality(reduced)
        self.assertEqual(quality["components"], 1)
        self.assertEqual(quality["boundary_edges"], 0)
        self.assertEqual(quality["nonmanifold_edges"], 0)
        small = vtk.vtkSphereSource()
        small.Update()
        _, unchanged = simplify_surface(small.GetOutput(), 150000)
        self.assertEqual(unchanged["input_triangles"], unchanged["output_triangles"])
        plane = vtk.vtkPlaneSource()
        plane.Update()
        with self.assertRaises(ValueError):
            simplify_surface(plane.GetOutput(), 150000)
