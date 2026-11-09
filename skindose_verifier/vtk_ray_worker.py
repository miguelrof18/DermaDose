"""Optional PythonSlicer worker: exact finite-line intersections on the full mesh."""
import json
from pathlib import Path
import sys
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray


def serve(folder):
    folder = Path(folder)
    vertices = np.load(folder/'vertices.npy', mmap_mode='r')
    faces = np.load(folder/'faces.npy', mmap_mode='r')
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(vertices, deep=False))
    cells = vtk.vtkCellArray()
    offsets = np.arange(len(faces)+1, dtype=np.int64)*3
    cells.SetData(numpy_to_vtkIdTypeArray(offsets, deep=False),
                  numpy_to_vtkIdTypeArray(faces.ravel(), deep=False))
    poly = vtk.vtkPolyData()
    poly.SetPoints(points)
    poly.SetPolys(cells)
    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(poly)
    locator.BuildLocator()
    print(json.dumps({'ready': True, 'vtk_version': vtk.vtkVersion.GetVTKVersion()}), flush=True)
    cell = vtk.vtkGenericCell()
    for line in sys.stdin:
        request = json.loads(line)
        if request.get('close'):
            return
        source = request['source']
        targets = np.load(folder/'targets.npy')
        result = np.full(len(targets), -1, dtype=np.int64)
        t, sub_id, cell_id = vtk.mutable(0.), vtk.mutable(0), vtk.mutable(0)
        hit, parametric = [0.,0.,0.], [0.,0.,0.]
        for index, target in enumerate(targets):
            vector = target - source
            distance = np.linalg.norm(vector)
            endpoint = target + vector * (1e-6/max(distance, 1e-12))
            if locator.IntersectWithLine(source, endpoint, 1e-7, t, hit, parametric, sub_id, cell_id, cell):
                if np.linalg.norm(np.asarray(hit)-source) <= distance+1e-7:
                    result[index] = int(cell_id)
        np.save(folder/'hits.npy', result)
        print(json.dumps({'done': True}), flush=True)


if __name__ == '__main__':
    serve(sys.argv[1])
