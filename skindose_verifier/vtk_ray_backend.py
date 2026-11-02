"""Optional full-resolution accelerator using an installed PythonSlicer.

No Slicer scene or registration is involved; only vertices, face indices and rays.
"""
import atexit
import hashlib
import json
import os
import re
from pathlib import Path
import subprocess
import tempfile
import threading
import numpy as np

_worker = None
_lock = threading.Lock()


def mesh_signature(mesh):
    return (id(mesh), hashlib.sha256(np.ascontiguousarray(mesh.vertices).data).digest(),
            hashlib.sha256(np.ascontiguousarray(mesh.faces).data).digest())


def find_python_slicer():
    configured = os.environ.get('DERMADOSE_PYTHONSLICER')
    if configured:
        path = Path(configured)
        return path if path.is_file() else None
    bases = [Path.home()/'AppData/Local/slicer.org']
    if os.environ.get('LOCALAPPDATA'):
        bases.append(Path(os.environ['LOCALAPPDATA'])/'slicer.org')
    candidates = sorted({p for base in bases for p in base.glob('3D Slicer */bin/PythonSlicer.exe')},
                        key=lambda p: tuple(int(n) for n in re.findall(r'\d+',p.parent.parent.name)), reverse=True)
    return candidates[0] if candidates else None


class VTKWorker:
    def __init__(self, mesh, executable):
        self.mesh = mesh
        self.signature = mesh_signature(mesh)
        self.temp = tempfile.TemporaryDirectory(prefix='dermadose_rays_', ignore_cleanup_errors=True)
        self.folder = Path(self.temp.name)
        self.process = None
        try:
            np.save(self.folder/'vertices.npy', np.asarray(mesh.vertices, dtype=np.float64))
            np.save(self.folder/'faces.npy', np.asarray(mesh.faces, dtype=np.int64))
            self.errors = (self.folder/'worker.log').open('w', encoding='utf-8')
            self.process = subprocess.Popen([str(executable), '-u', str(Path(__file__).with_name('vtk_ray_worker.py')), str(self.folder)],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.errors, text=True, encoding='utf-8',
                creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
            self.vtk_version = self.response('ready').get('vtk_version')
        except Exception:
            self.close()
            raise

    def response(self, key):
        line = self.process.stdout.readline()
        result = json.loads(line) if line else {}
        if not result.get(key):
            raise RuntimeError('El acelerador VTK no devolvió una respuesta válida.')
        return result

    def intersect(self, source, targets):
        np.save(self.folder/'targets.npy', targets)
        self.process.stdin.write(json.dumps({'source': np.asarray(source).tolist()})+'\n')
        self.process.stdin.flush()
        self.response('done')
        return np.load(self.folder/'hits.npy')

    def close(self):
        if self.process is not None:
            if self.process.poll() is None:
                try:
                    self.process.stdin.write('{"close":true}\n')
                    self.process.stdin.flush()
                except (OSError, ValueError):
                    pass
            if self.process.stdin:
                self.process.stdin.close()
            try:
                self.process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait(timeout=10)
            if self.process.stdout:
                self.process.stdout.close()
        if hasattr(self, 'errors'):
            self.errors.close()
        # TemporaryDirectory owns only files created above, under the temp root.
        if self.folder.resolve().parent == Path(tempfile.gettempdir()).resolve():
            self.temp.cleanup()


def close_worker():
    global _worker
    if _worker is not None:
        _worker.close()
        _worker = None


atexit.register(close_worker)


def intersect(mesh, source, targets):
    global _worker
    executable = find_python_slicer()
    if executable is None:
        return None
    with _lock:
        if _worker is None or _worker.signature != mesh_signature(mesh):
            close_worker()
            _worker = VTKWorker(mesh, executable)
        return _worker.intersect(source, targets), _worker.vtk_version
