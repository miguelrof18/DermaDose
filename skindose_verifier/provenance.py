"""Environment and source fingerprints for reproducible research outputs."""
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path
import platform


def runtime_provenance() -> dict:
    root = Path(__file__).parent
    paths = sorted(root.glob("*.py")) + [root/"data/visor_3d_template.html"]
    return {"python": platform.python_version(), "platform": platform.platform(),
            "dependencies": {name: version(name) for name in ("numpy", "scipy", "trimesh", "rtree", "pydicom", "plotly", "matplotlib")},
            "source_sha256": {str(path.relative_to(root)): sha256(path.read_bytes()).hexdigest() for path in paths}}
