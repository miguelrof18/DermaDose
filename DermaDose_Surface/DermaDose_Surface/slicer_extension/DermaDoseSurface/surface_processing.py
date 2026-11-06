"""Standalone VTK validation, shared by Slicer and its extraction worker."""
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk


def triangle_statistics(poly):
    """Dimensionless mean-ratio quality: 1 equilateral, 0 degenerate."""
    points = vtk_to_numpy(poly.GetPoints().GetData())
    faces = vtk_to_numpy(poly.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    qualities, edges, areas = [], [], []
    for start in range(0, len(faces), 100000):
        p = points[faces[start:start+100000]].astype(np.float64)
        e = np.stack((p[:,1]-p[:,0], p[:,2]-p[:,1], p[:,0]-p[:,2]), axis=1)
        squared = np.einsum('ijk,ijk->ij', e, e)
        area = .5*np.linalg.norm(np.cross(e[:,0], -e[:,2]), axis=1)
        quality = np.divide(4*np.sqrt(3)*area, squared.sum(axis=1),
                            out=np.zeros(len(area)), where=squared.sum(axis=1)>0)
        qualities.append(quality)
        edges.append(np.sqrt(squared).ravel())
        areas.append(area)
    quality, lengths, area = np.concatenate(qualities), np.concatenate(edges), np.concatenate(areas)
    return dict(mean_ratio_p01=float(np.percentile(quality,1)),
                mean_ratio_p05=float(np.percentile(quality,5)),
                mean_ratio_median=float(np.median(quality)),
                slender_face_fraction=float(np.mean(quality < .1)),
                edge_length_p95_mm=float(np.percentile(lengths,95)),
                edge_length_max_mm=float(lengths.max()), area_mm2=float(area.sum()))

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
            # A local origin avoids cancellation for large DICOM offsets.
            volume += float(np.einsum("ij,ij->", vertices[:, 0] - points[0], cross) / 6.0)
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
    """Bidirectional vertices and face centers; not a continuous Hausdorff bound."""
    distance = vtk.vtkDistancePolyDataFilter()
    distance.SetInputData(0, source)
    distance.SetInputData(1, candidate)
    distance.SignedDistanceOff()
    distance.ComputeSecondDistanceOn()
    distance.ComputeCellCenterDistanceOn()
    distance.Update()
    values = []
    for i in (0, 1):
        output = distance.GetOutput(i)
        # VTK evaluates both vertices and cell centers with the same locator.
        # Read both arrays; computing a second center-only filter is unnecessary.
        for attributes in (output.GetPointData(), output.GetCellData()):
            array = attributes.GetArray('Distance')
            if array is None:
                raise ValueError('No se pudo comprobar la desviación de la superficie.')
            v = np.abs(vtk_to_numpy(array))
            if not len(v) or not np.isfinite(v).all():
                raise ValueError('Distancia de superficie no finita.')
            values.append(v)
    return float(max(v.max() for v in values))


def simplify_surface(poly, target_triangles=150000, progress=None,
                     max_deviation_mm=1.0, max_volume_change=0.01,
                     minimum_triangles=None, maximum_triangles=None, optimization_mode="fast"):
    """Accept a count band only when topology and geometric tolerances pass.

    No smoothing, hole filling, subdivision or coordinate-system changes.
    The unchanged source remains available for inspection if reduction fails.
    """
    if optimization_mode not in ("fast", "compact"):
        raise ValueError("Prioridad de optimización desconocida.")
    target_triangles = int(target_triangles)
    if target_triangles < 100:
        raise ValueError("El objetivo mínimo es de 100 triángulos.")
    if not np.isfinite(max_deviation_mm) or max_deviation_mm <= 0:
        raise ValueError("La tolerancia geométrica debe ser positiva y finita.")
    lower = int(minimum_triangles if minimum_triangles is not None else round(target_triangles * 2 / 3))
    upper = int(maximum_triangles if maximum_triangles is not None else round(target_triangles * 4 / 3))
    if not (100 <= lower <= target_triangles <= upper):
        raise ValueError("Margen de caras inválido.")
    if not np.isfinite(max_volume_change) or not 0 <= max_volume_change < 1:
        raise ValueError("Tolerancia de volumen inválida.")
    source = clean_surface(poly)
    def orient(candidate):
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputData(candidate)
        normals.SplittingOff()
        normals.ConsistencyOn()
        normals.AutoOrientNormalsOn()
        normals.Update()
        return clean_surface(normals.GetOutput())
    before = surface_quality(source)
    if not valid_surface(before):
        raise ValueError("Superficie original inválida: " + str(before)
                         + ". Corrija la segmentación antes de simplificar.")
    source = orient(source)
    before = surface_quality(source)
    if not valid_surface(before):
        raise ValueError("No se pudo orientar la superficie original.")
    count = before["triangles"]
    source_stats = triangle_statistics(source)
    euler_before = source.GetNumberOfPoints() - count // 2
    result, chosen_method, attempts = source, "unchanged", []
    metrics = dict(max_vertex_surface_distance_mm=0.0, relative_volume_change=0.0,
                   relative_area_change=0.0)
    accepted = count <= upper
    # Existing meshes in the band are already usable; don't move anatomy merely
    # to shave off additional faces. Extraction starts from the full-resolution mesh.
    best_score = None
    if not accepted:
        budgets = list(dict.fromkeys((target_triangles, upper) if optimization_mode == "fast"
                                     else (target_triangles, lower, upper)))
        for budget in budgets:
            # A coarser accepted mesh makes a finer budget unnecessary.
            if best_score is not None and budget > result.GetNumberOfPolys():
                continue
            methods = (("topology_preserving", "quadric_normals_regularized") if optimization_mode == "fast"
                       else ("quadric_normals_regularized", "topology_preserving"))
            for method in methods:
                if progress:
                    progress(f"Optimizando: {method}, {budget} caras; margen {lower}–{upper}")
                settings = {}
                if method == "topology_preserving":
                    decimator = vtk.vtkDecimatePro()
                    decimator.PreserveTopologyOn()
                    decimator.SplittingOff()
                    decimator.BoundaryVertexDeletionOff()
                    decimator.SetFeatureAngle(180.0)
                else:
                    decimator = vtk.vtkQuadricDecimation()
                    decimator.VolumePreservationOn()
                    decimator.AttributeErrorMetricOn()
                    decimator.ScalarsAttributeOff()
                    decimator.VectorsAttributeOff()
                    decimator.TCoordsAttributeOff()
                    decimator.TensorsAttributeOff()
                    decimator.NormalsAttributeOn()
                    normal_weight = .5*np.sqrt(4*source_stats["area_mm2"]/(np.sqrt(3)*budget))
                    decimator.SetNormalsWeight(float(normal_weight))
                    settings["normal_weight"] = float(normal_weight)
                    settings["regularization_enabled"] = hasattr(decimator,"RegularizeOn")
                    if settings["regularization_enabled"]:
                        decimator.RegularizeOn()
                        # VTK's probabilistic quadric regularization, not smoothing.
                        regularization = min(.05, max_deviation_mm*.25)
                        decimator.SetRegularization(regularization)
                        settings["regularization"] = regularization
                decimator.SetInputData(source)
                decimator.SetTargetReduction(1.0 - budget / count)
                decimator.Update()
                if progress:
                    progress(f"Verificando topología de {method} ({budget:,} caras)")
                candidate = orient(clean_surface(decimator.GetOutput()))
                if candidate.GetNumberOfPoints() == 0:
                    attempts.append(dict(method=method, requested_triangles=budget,
                                         accepted=False, rejection="empty_output"))
                    continue
                candidate.GetPoints().SetData(numpy_to_vtk(
                    vtk_to_numpy(candidate.GetPoints().GetData()).astype(np.float32), deep=True))
                candidate = clean_surface(candidate)
                quality = surface_quality(candidate)
                attempt = dict(method=method, requested_triangles=budget, settings=settings, **quality)
                topology_ok = (valid_surface(quality)
                    and candidate.GetNumberOfPoints() - quality["triangles"] // 2 == euler_before)
                attempt["topology_preserved"] = topology_ok
                attempt["within_count_range"] = lower <= quality["triangles"] <= upper
                candidate_ok = False
                if topology_ok and attempt["within_count_range"]:
                    stats = triangle_statistics(candidate)
                    volume_change = abs(abs(quality["signed_volume_mm3"]) /
                                        abs(before["signed_volume_mm3"]) - 1)
                    area_change = abs(stats["area_mm2"]/source_stats["area_mm2"]-1)
                    shape_ok = stats["mean_ratio_p01"] >= min(.05, source_stats["mean_ratio_p01"])
                    if volume_change <= max_volume_change and area_change <= .05 and shape_ok:
                        if progress:
                            progress(f"Comprobando desviación geométrica ({budget:,} caras)")
                        deviation = shape_difference(source, candidate)
                        candidate_ok = deviation <= max_deviation_mm
                    else:
                        deviation = None
                    attempt.update(max_vertex_surface_distance_mm=deviation,
                                   relative_volume_change=volume_change, relative_area_change=area_change,
                                   triangle_quality=stats, triangle_shape_passed=bool(shape_ok))
                    if candidate_ok:
                        # Counts within 500 faces are treated equally: prefer fewer
                        # sliver triangles, then smaller geometric error.
                        score = (int(np.ceil(quality["triangles"]/500)),
                                 stats["slender_face_fraction"], deviation)
                        if best_score is None or score < best_score:
                            best_score = score
                            result, chosen_method = candidate, method
                            metrics = dict(max_vertex_surface_distance_mm=deviation,
                                           relative_volume_change=volume_change, relative_area_change=area_change)
                        accepted = True
                attempt["accepted"] = bool(candidate_ok)
                attempts.append(attempt)
                if progress and not candidate_ok:
                    progress(f"Candidato descartado: {quality['triangles']:,} caras; "
                             f"topología={topology_ok}; desviación={attempt.get('max_vertex_surface_distance_mm')}; "
                             f"cambio de área={attempt.get('relative_area_change')}")
                if candidate_ok and optimization_mode == "fast":
                    break
            if accepted and optimization_mode == "fast":
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
        source_triangle_quality=source_stats,
        triangle_quality=source_stats if result is source else triangle_statistics(result),
        surface_quality=after, optimization_mode=optimization_mode,
        selection_policy=("first validated candidate; alternatives only on failure" if optimization_mode == "fast"
                          else "fewest compliant faces among searched candidates; 500-face buckets then shape and distance"),
        globally_optimal=False, max_deviation_tolerance_mm=max_deviation_mm,
        relative_volume_tolerance=max_volume_change, relative_area_tolerance=.05,
        distance_measure="bidirectional vertices and all face centroids to opposing surface; sampled, not continuous Hausdorff",
        self_intersections_checked=False, precision_validated=False, **metrics,
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
    if quality["triangles"] != poly.GetNumberOfPolys():
        raise ValueError("El STL guardado cambia el número de caras.")
    # STL preserves triangle order. A direct correspondence provides an upper
    # bound for every point inside each triangle, without building two BVHs.
    incoming = reader.GetOutput()
    points_a = vtk_to_numpy(poly.GetPoints().GetData())
    points_b = vtk_to_numpy(incoming.GetPoints().GetData())
    faces_a = vtk_to_numpy(poly.GetPolys().GetConnectivityArray()).reshape(-1,3)
    faces_b = vtk_to_numpy(incoming.GetPolys().GetConnectivityArray()).reshape(-1,3)
    roundtrip_error = 0.0
    for start in range(0, len(faces_a), 100000):
        a = points_a[faces_a[start:start+100000]].astype(np.float64)
        b = points_b[faces_b[start:start+100000]].astype(np.float64)
        roundtrip_error = max(roundtrip_error, float(np.linalg.norm(a-b,axis=2).max()))
    check_method = "corresponding_triangle_vertices"
    if roundtrip_error > .001:
        # Readers may reorder triangles. Retain the general geometric fallback.
        roundtrip_error = shape_difference(poly, incoming)
        check_method = "bidirectional_vertices_and_centroids"
    if roundtrip_error > .001:
        raise ValueError("La precisión STL altera la superficie más de 0,001 mm.")
    quality["stl_roundtrip_sampled_distance_mm"] = roundtrip_error
    quality["stl_roundtrip_check_method"] = check_method
    return quality



