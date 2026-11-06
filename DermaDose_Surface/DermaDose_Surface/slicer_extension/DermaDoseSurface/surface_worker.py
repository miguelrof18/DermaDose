"""Standalone extraction; no Qt or MRML objects are accessed here."""
import json
from pathlib import Path
import sys
import time
import traceback
import numpy as np
from scipy import ndimage
import vtk
from vtk.util.numpy_support import numpy_to_vtk
from surface_processing import simplify_surface


def largest_body(values, threshold, progress=print):
    if values.ndim != 3 or not np.isfinite(threshold):
        raise ValueError('Se requiere un TC escalar 3D y un umbral finito.')
    if values.size > 160_000_000:
        raise ValueError('TC de más de 160 millones de vóxeles. Recorte una ROI o remuestree el TC antes de extraer.')
    progress('Aplicando umbral')
    mask = values >= threshold
    if not np.any(mask):
        raise ValueError('El umbral no selecciona ningún vóxel.')
    progress('Buscando el componente corporal principal')
    labels, count = ndimage.label(mask)
    del mask
    sizes = np.zeros(count + 1, dtype=np.int64)
    flat = labels.reshape(-1)
    for start in range(0, len(flat), 1_000_000):
        sizes += np.bincount(flat[start:start+1_000_000], minlength=count + 1)
    sizes[0] = 0
    body = labels == int(np.argmax(sizes))
    del flat, labels, sizes
    progress('Rellenando cavidades internas')
    return ndimage.binary_fill_holes(body).astype(np.uint8)


def make_surface(body, ijk_to_ras):
    image = vtk.vtkImageData()
    image.SetDimensions(*body.shape[::-1])
    image.GetPointData().SetScalars(numpy_to_vtk(body.ravel(), deep=False))
    pad = vtk.vtkImageConstantPad()
    pad.SetInputData(image)
    nx, ny, nz = body.shape[::-1]
    pad.SetOutputWholeExtent(-1, nx, -1, ny, -1, nz)
    pad.SetConstant(0)
    contour = vtk.vtkFlyingEdges3D()
    contour.SetInputConnection(pad.GetOutputPort())
    contour.SetValue(0, 0.5)
    contour.ComputeNormalsOff()
    contour.Update()
    matrix = vtk.vtkMatrix4x4()
    for i in range(4):
        for j in range(4):
            matrix.SetElement(i, j, ijk_to_ras[i][j])
    transform = vtk.vtkTransform()
    transform.SetMatrix(matrix)
    apply = vtk.vtkTransformPolyDataFilter()
    apply.SetInputConnection(contour.GetOutputPort())
    apply.SetTransform(transform)
    apply.Update()
    poly = vtk.vtkPolyData()
    poly.DeepCopy(apply.GetOutput())
    return poly


def run(folder):
    folder = Path(folder)
    job = json.loads((folder/'job.json').read_text(encoding='utf-8'))
    def progress(message):
        print(message, flush=True)
    start = time.perf_counter()
    if job.get('mode') == 'review':
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(str(folder/'original_surface.vtp'))
        reader.Update()
        poly = reader.GetOutput()
    else:
        values = np.load(folder/'input.npy', mmap_mode='r')
        body = largest_body(values, job['threshold_hu'], progress)
        del values
        progress('Construyendo superficie cerrada en coordenadas RAS')
        poly = make_surface(body, job['ijk_to_ras'])
        original_writer = vtk.vtkXMLPolyDataWriter()
        original_writer.SetCompressorTypeToNone()
        original_writer.EncodeAppendedDataOff()
        original_writer.SetFileName(str(folder/'original_surface.vtp'))
        original_writer.SetInputData(poly)
        if not original_writer.Write() or original_writer.GetErrorCode():
            raise OSError('No se pudo conservar la superficie original.')
        np.save(folder/'body.npy', body)
    poly, decimation = simplify_surface(poly, job['target_triangles'], progress,
                                       max_deviation_mm=job.get('max_deviation_mm', 1.0),
                                       optimization_mode=job.get('optimization_mode','fast'),
                                       fallback_step_triangles=100000)
    progress('Guardando el resultado para Slicer')
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetCompressorTypeToNone()
    writer.EncodeAppendedDataOff()
    writer.SetFileName(str(folder/'surface.vtp'))
    writer.SetInputData(poly)
    if not writer.Write() or writer.GetErrorCode():
        raise OSError('No se pudo guardar la superficie temporal.')
    decimation['elapsed_seconds'] = time.perf_counter() - start
    (folder/'result.json').write_text(json.dumps(decimation), encoding='utf-8')
    progress(f"Extracción terminada en {decimation['elapsed_seconds']:.1f} s")


if __name__ == '__main__':
    try:
        run(sys.argv[1])
    except Exception:
        traceback.print_exc()
        sys.exit(1)

