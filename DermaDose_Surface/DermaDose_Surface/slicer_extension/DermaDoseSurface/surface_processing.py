"""Standalone VTK validation, shared by Slicer and its extraction worker."""
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

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


def simplify_surface(poly, target_triangles=150000, progress=None):
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
        budgets = [int(target_triangles)]
    chosen_method = "unchanged"
    found = False
    for budget in budgets:
        if budget >= count:
            continue
        # Try both available VTK decimators; accept only independently checked outputs.
        for method in ("vtkDecimatePro", "vtkQuadricDecimation"):
            if progress:
                progress(f'Simplificando: {method}, {count} a {budget} triángulos')
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


