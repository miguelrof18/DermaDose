"""Generate independent comparison meshes from the same original surface."""
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
import vtk
from surface_processing import simplify_surface, write_stl_checked

def export_pair(source_path, destination):
    source_path, destination = Path(source_path), Path(destination)
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(source_path))
    reader.Update()
    source = reader.GetOutput()
    if source.GetNumberOfPolys() < 300000:
        raise ValueError("El original no tiene 300 000 caras. No se añadirán caras artificiales.")
    names = ["surface_150k.stl", "surface_150k.json",
             "surface_300k.stl", "surface_300k.json", "comparison.json"]
    if any((destination / name).exists() for name in names):
        raise FileExistsError("Ya existe un archivo de comparación. Elija otra carpeta.")
    source_hash = hashlib.sha256(source_path.read_bytes()).hexdigest()
    created = []
    with tempfile.TemporaryDirectory(dir=str(destination)) as temporary:
        temporary = Path(temporary)
        results = {}
        for target, minimum, maximum in ((150000,100000,200000),(300000,270000,330000)):
            poly, metadata = simplify_surface(source, target, print,
                minimum_triangles=minimum, maximum_triangles=maximum)
            if not metadata["export_allowed"] or not minimum <= metadata["output_triangles"] <= maximum:
                raise ValueError(f"No se obtuvo una malla válida para {target}: {metadata}")
            name = f"surface_{target//1000}k"
            metadata["stl_quality"] = write_stl_checked(poly, temporary/(name+".stl"))
            metadata.update(coordinate_system="RAS", length_unit="mm", source_sha256=source_hash,
                mesh_sha256=hashlib.sha256((temporary/(name+".stl")).read_bytes()).hexdigest())
            (temporary/(name+".json")).write_text(json.dumps(metadata,indent=2),encoding="utf-8")
            results[name] = metadata
        manifest = dict(source_sha256=source_hash, meshes=results, psd_compared=False,
            instructions="Calcular ambas con el mismo RDSR, configuración y registro. "
            "Comparar PSD, diferencia relativa y posición del máximo. "
            "300 000 caras tampoco constituyen una referencia exacta.")
        (temporary/"comparison.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
        try:
            for name in names:
                with (destination/name).open("xb") as output:
                    created.append(destination/name)
                    with (temporary/name).open("rb") as incoming:
                        shutil.copyfileobj(incoming, output)
        except Exception:
            for path in created:
                path.unlink(missing_ok=True)
            raise
    return manifest

if __name__ == "__main__":
    export_pair(sys.argv[1], sys.argv[2])
