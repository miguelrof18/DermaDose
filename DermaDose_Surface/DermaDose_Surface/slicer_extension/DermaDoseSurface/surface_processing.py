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


def shape_difference(source, candidate):
    """Bidirectional vertex-to-surface distances, not a continuous Hausdorff bound."""
    distance = vtk.vtkDistancePolyDataFilter()
    distance.SetInputData(0, source)
    distance.SetInputData(1, candidate)
    distance.SignedDistanceOff()
    distance.ComputeSecondDistanceOn()
    distance.Update()
    values = [np.abs(vtk_to_numpy(distance.GetOutput(i).GetPointData().GetArray("Distance")))
              for i in (0, 1)]
    if any(not len(v) or not np.isfinite(v).all() for v in values):
        raise ValueError("No se pudo comprobar la desviación de la superficie.")
    return float(max(v.max() for v in values))


def simplify_surface(poly, target_triangles=150000, progress=None,
                     max_deviation_mm=1.0, max_volume_change=0.01,
                     minimum_triangles=None, maximum_triangles=None):
    """Accept a count band only when topology and geometric tolerances pass.

    No smoothing, hole filling, subdivision or coordinate-system changes.
    The unchanged source remains available for inspection if reduction fails.
    """
    target_triangles = int(target_triangles)
    if target_triangles < 100:
        raise ValueError("El objetivo mínimo es de 100 triángulos.")
    if not np.isfinite(max_deviation_mm) or max_deviation_mm <= 0:
        raise ValueError("La tolerancia geométrica debe ser positiva y finita.")
    lower = int(minimum_triangles if minimum_triangles is not None else round(target_triangles * 2 / 3))
    upper = int(maximum_triangles if maximum_triangles is not None else round(target_triangles * 4 / 3))
    source = clean_surface(poly)
    def orient(candidate):
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputData(candidate)
        normals.SplittingOff()
        normals.ConsistencyOn()
        normals.AutoOrientNormalsOn()
        normals.Update()
        return clean_surface(normals.GetOutput())
    source = orient(source)
    before = surface_quality(source)
    if not valid_surface(before):
        raise ValueError("Superficie original inválida: " + str(before)
                         + ". Corrija la segmentación antes de simplificar.")
    count = before["triangles"]
    euler_before = source.GetNumberOfPoints() - count // 2
    result, chosen_method, attempts = source, "unchanged", []
    metrics = dict(max_vertex_surface_distance_mm=0.0, relative_volume_change=0.0)
    accepted = count <= upper
    if not accepted:
        # Raising the feature angle releases voxel-scale facets without allowing
        # topology changes. Every candidate is independently checked afterwards.
        for method, budget in (("vtkDecimatePro", target_triangles),
                               ("vtkQuadricDecimation", target_triangles),
                               ("vtkDecimatePro", upper),
                               ("vtkQuadricDecimation", upper)):
            if progress:
                progress(f"Simplificando: {method}, objetivo {budget}; margen {lower}–{upper}")
            decimator = getattr(vtk, method)()
            decimator.SetInputData(source)
            decimator.SetTargetReduction(1.0 - budget / count)
            if method == "vtkDecimatePro":
                decimator.PreserveTopologyOn()
                decimator.SplittingOff()
                decimator.BoundaryVertexDeletionOff()
                decimator.SetFeatureAngle(180.0)
            elif hasattr(decimator, "VolumePreservationOn"):
                decimator.VolumePreservationOn()
            decimator.Update()
            candidate = orient(clean_surface(decimator.GetOutput()))
            # STL stores float32 coordinates. Validate that precision before acceptance.
            candidate.GetPoints().SetData(
                __import__("vtk.util.numpy_support", fromlist=["numpy_to_vtk"]).numpy_to_vtk(
                    vtk_to_numpy(candidate.GetPoints().GetData()).astype(np.float32), deep=True))
            candidate = clean_surface(candidate)
            quality = surface_quality(candidate)
            attempt = dict(method=method, requested_triangles=budget, **quality)
            topology_ok = (valid_surface(quality)
                and candidate.GetNumberOfPoints() - quality["triangles"] // 2 == euler_before)
            attempt["topology_preserved"] = topology_ok
            attempt["within_count_range"] = lower <= quality["triangles"] <= upper
            if topology_ok and attempt["within_count_range"]:
                deviation = shape_difference(source, candidate)
                volume_change = abs(abs(quality["signed_volume_mm3"]) /
                                    abs(before["signed_volume_mm3"]) - 1)
                attempt.update(max_vertex_surface_distance_mm=deviation,
                               relative_volume_change=volume_change)
                accepted = deviation <= max_deviation_mm and volume_change <= max_volume_change
                if accepted:
                    result, chosen_method = candidate, method
                    metrics = dict(max_vertex_surface_distance_mm=deviation,
                                   relative_volume_change=volume_change)
            attempt["accepted"] = bool(accepted)
            attempts.append(attempt)
            if accepted:
                break
    after = surface_quality(result)
    return result, dict(
        method=chosen_method, target_triangles=target_triangles,
        minimum_triangles=lower, maximum_triangles=upper,
        input_triangles=count, output_triangles=after["triangles"],
        within_range=lower <= after["triangles"] <= upper,
        target_reached=bool(accepted), export_allowed=bool(accepted),
        below_range_without_subdivision=count < lower,
        fallback_original=not accepted, attempts=attempts,
        max_deviation_tolerance_mm=max_deviation_mm,
        relative_volume_tolerance=max_volume_change,
        distance_measure="bidirectional vertices to opposing surface; not continuous Hausdorff",
        precision_validated=False, **metrics,
    )


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


