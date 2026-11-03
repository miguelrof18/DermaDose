"""3D Slicer module: CT external contour for research use in DermaDose."""

import hashlib
import json
import tempfile
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


class DermaDoseSurface(ScriptedLoadableModule):
    def __init__(self, parent):
        super().__init__(parent)
        parent.title = "DermaDose Surface"
        parent.categories = ["DermaDose"]
        parent.contributors = ["DermaDose"]
        parent.helpText = (
            "Extract an external CT contour, inspect it, and export a RAS/mm STL "
            "with an audit JSON. This is not a skin-layer segmentation."
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

        self.threshold = qt.QDoubleSpinBox()
        self.threshold.minimum = -1500
        self.threshold.maximum = 3000
        self.threshold.value = -350
        self.threshold.suffix = " HU"
        self.threshold.toolTip = "Umbral inicial orientativo; revise el contorno en los cortes."
        form.addRow("Umbral mínimo:", self.threshold)

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

    def onExtract(self):
        try:
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
        slicer.util.selectModule("SegmentEditor")
        editor = slicer.modules.segmenteditor.widgetRepresentation().self().editor
        editor.setSegmentationNode(self.logic.segmentationNode)
        editor.setSourceVolumeNode(self.logic.volumeNode)
        editor.setCurrentSegmentID(self.logic.segmentId)
        self.acknowledge.checked = False

    def onCheck(self):
        try:
            result = self.logic.inspect()
            lines = [
                f"Triángulos: {result['triangles']}",
                f"Componentes de superficie: {result['components']}",
                f"Aristas abiertas: {result['boundary_edges']}",
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
        if self.segmentationNode:
            slicer.mrmlScene.RemoveNode(self.segmentationNode)
        seg = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", "DermaDose external contour")
        try:
            seg.CreateDefaultDisplayNodes()
            seg.SetReferenceImageGeometryParameterFromVolumeNode(volume)
            segment_id = seg.GetSegmentation().AddEmptySegment("ExternalContour", "External contour")
            slicer.util.updateSegmentBinaryLabelmapFromArray(body, seg, segment_id, volume)
            if not seg.CreateClosedSurfaceRepresentation():
                raise RuntimeError("No se pudo generar la superficie cerrada.")
            seg.GetDisplayNode().SetAllSegmentsVisibility(True)
            self.segmentationNode, self.volumeNode = seg, volume
            self.segmentId, self.thresholdHU = segment_id, float(threshold_hu)
            slicer.util.resetThreeDViews()
            return seg
        except Exception:
            slicer.mrmlScene.RemoveNode(seg)
            raise

    def _surface(self):
        if not self.segmentationNode or not self.segmentId:
            raise ValueError("Extraiga primero la superficie del TC.")
        if self.volumeNode.GetParentTransformNode() or self.segmentationNode.GetParentTransformNode():
            raise ValueError("Endurezca las transformaciones antes de exportar.")
        self.segmentationNode.CreateClosedSurfaceRepresentation()
        poly = vtk.vtkPolyData()
        if not self.segmentationNode.GetClosedSurfaceRepresentation(self.segmentId, poly):
            raise ValueError("No se pudo obtener la superficie del segmento.")
        if poly.GetNumberOfPolys() == 0:
            raise ValueError("El contorno no contiene triángulos.")
        return poly

    def inspect(self):
        poly = self._surface()
        connectivity = vtk.vtkPolyDataConnectivityFilter()
        connectivity.SetInputData(poly)
        connectivity.SetExtractionModeToAllRegions()
        connectivity.Update()
        edges = vtk.vtkFeatureEdges()
        edges.SetInputData(poly)
        edges.BoundaryEdgesOn()
        edges.FeatureEdgesOff()
        edges.ManifoldEdgesOff()
        edges.NonManifoldEdgesOff()
        edges.Update()
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
        if connectivity.GetNumberOfExtractedRegions() != 1:
            warnings.append("Hay varios componentes de superficie; revise objetos extraños.")
        if edges.GetOutput().GetNumberOfCells():
            warnings.append("Hay aristas abiertas; revise la integridad de la malla.")
        return {
            "triangles": int(poly.GetNumberOfPolys()),
            "components": int(connectivity.GetNumberOfExtractedRegions()),
            "boundary_edges": int(edges.GetOutput().GetNumberOfCells()),
            "warnings": warnings,
        }

    def export(self, path):
        path = Path(path)
        if path.suffix.lower() != ".stl":
            raise ValueError("Seleccione un nombre terminado en .stl.")
        if path.exists() or path.with_suffix(".json").exists():
            raise FileExistsError("Ya existe el STL o su JSON; elija un nombre nuevo.")
        qa = self.inspect()
        poly = self._surface()
        triangles = vtk.vtkTriangleFilter()
        triangles.SetInputData(poly)
        triangles.Update()
        writer = vtk.vtkSTLWriter()
        writer.SetFileName(str(path))
        writer.SetFileTypeToBinary()
        writer.SetHeader("DermaDose SPACE=RAS; coordinates in mm")
        writer.SetInputConnection(triangles.GetOutputPort())
        if not writer.Write():
            raise OSError("No se pudo escribir el STL.")
        matrix = vtk.vtkMatrix4x4()
        self.volumeNode.GetIJKToRASMatrix(matrix)
        metadata = {
            "schema_version": 1,
            "surface_definition": "external CT contour, not histological skin",
            "mesh_file": path.name,
            "mesh_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "coordinate_system": "RAS",
            "length_unit": "mm",
            "ijk_to_ras_mm": [[matrix.GetElement(i, j) for j in range(4)] for i in range(4)],
            "threshold_minimum_hu": self.thresholdHU,
            "processing": ["threshold", "largest_6_connected_component", "fill_enclosed_holes"],
            "manual_edits": "Review the scene/segmentation; manual changes are not automatically recorded",
            "quality": qa,
            "rdsr_registration": "Not established; verify geometry separately before dose calculation",
        }
        try:
            path.with_suffix(".json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        except Exception:
            path.unlink(missing_ok=True)
            raise
        return metadata


class DermaDoseSurfaceTest(ScriptedLoadableModuleTest):
    def runTest(self):
        self.setUp()
        self.testSyntheticContour()

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
