"""Central-projection visibility against brute-force Möller–Trumbore ray casting.

Pure NumPy: closed, outward-oriented meshes with strong self-occlusion.
"""
import numpy as np
import pytest

from skindose_verifier.raytrace import visible_faces


def icosphere(levels, radius):
    t = (1+5**.5)/2
    v = np.array([[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t],
                  [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]], float)
    f = np.array([[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4],
                  [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8],
                  [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]])
    v /= np.linalg.norm(v, axis=1)[:, None]
    for _ in range(levels):
        vertices, cache, faces = list(v), {}, []

        def middle(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                m = (vertices[a]+vertices[b])/2
                vertices.append(m/np.linalg.norm(m))
                cache[key] = len(vertices)-1
            return cache[key]
        for a, b, c in f:
            ab, bc, ca = middle(a, b), middle(b, c), middle(c, a)
            faces += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
        v, f = np.array(vertices), np.array(faces)
    return v*radius, f


def torus(major=100., minor=40., nu=60, nv=30):
    u, w = np.meshgrid(np.linspace(0, 2*np.pi, nu, endpoint=False),
                       np.linspace(0, 2*np.pi, nv, endpoint=False), indexing="ij")
    v = np.column_stack(((major+minor*np.cos(w)).ravel()*np.cos(u).ravel(),
                         (major+minor*np.cos(w)).ravel()*np.sin(u).ravel(), (minor*np.sin(w)).ravel()))
    index = lambda i, j: (i % nu)*nv + j % nv
    f = [[index(i, j), index(i+1, j), index(i+1, j+1)] for i in range(nu) for j in range(nv)]
    f += [[index(i, j), index(i+1, j+1), index(i, j+1)] for i in range(nu) for j in range(nv)]
    return v, np.array(f)


def oriented(v, f):
    volume = np.einsum("ij,ij->i", v[f[:, 0]], np.cross(v[f[:, 1]], v[f[:, 2]])).sum()/6
    f = f if volume > 0 else f[:, ::-1]
    n = np.cross(v[f[:, 1]]-v[f[:, 0]], v[f[:, 2]]-v[f[:, 0]])
    return f, n/np.linalg.norm(n, axis=1)[:, None]


def brute_force(v, f, source, targets):
    tri = v[f]
    e1, e2 = tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]
    out = np.zeros(len(targets), bool)
    for k, (face, point) in enumerate(zip(targets, tri[targets].mean(axis=1))):
        d = point-source
        h = np.cross(d, e2)
        a = np.einsum("ij,ij->i", e1, h)
        ok = np.abs(a) > 1e-14
        inv = np.where(ok, 1/np.where(ok, a, 1), 0)
        s = source-tri[:, 0]
        u = inv*np.einsum("ij,ij->i", s, h)
        q = np.cross(s, e1)
        w = inv*(q @ d)
        t = inv*np.einsum("ij,ij->i", e2, q)
        hit = ok & (u >= -1e-12) & (w >= -1e-12) & (u+w <= 1+1e-12) & (t > 1e-9)
        t = np.where(hit, t, np.inf)
        out[k] = hit.any() and (np.argmin(t) == face or abs(t[face]-t.min()) <= 1e-7)
    return out


def scenes():
    yield "sphere", *icosphere(3, 100.)
    a, fa = icosphere(3, 60.)
    b, fb = icosphere(2, 40.)
    yield "two_spheres", np.vstack((a, b+[30., 20., 150.])), np.vstack((fa, fb+len(a)))
    v, f = icosphere(3, 100.)
    yield "bumpy", v*(1+.15*np.sin(5*v[:, [0]]/100)*np.cos(4*v[:, [1]]/100)), f
    yield "torus", *torus()


@pytest.mark.parametrize("name,vertices,faces", list(scenes()), ids=lambda x: x if isinstance(x, str) else "")
def test_projective_visibility_matches_brute_force(name, vertices, faces):
    rng = np.random.default_rng(7)
    faces, normals = oriented(vertices, faces)
    for _ in range(4):
        source = rng.normal(size=3)
        source *= rng.uniform(400, 900)/np.linalg.norm(source)
        rays = vertices[faces].mean(axis=1)-source
        front = np.flatnonzero(np.einsum("ij,ij->i", normals, -rays) > 1e-9)
        axis = -source/np.linalg.norm(source)
        cosine = rays[front] @ axis/np.linalg.norm(rays[front], axis=1)
        targets = front[cosine > np.cos(np.radians(rng.uniform(8, 25)))]
        result = visible_faces(vertices, faces, normals, source, targets)
        assert result is not None
        visible, stats = result
        np.testing.assert_array_equal(visible, brute_force(vertices, faces, source, targets))
        assert stats["rays_traced"] == len(targets)
        assert stats["backend"] == "numpy_central_projection_cpu"


def test_occluded_sphere_is_hidden_and_focal_plane_falls_back():
    v, f = icosphere(2, 20.)
    near, far = v+[0., 0., 100.], v
    vertices = np.vstack((near, far))
    faces, normals = oriented(vertices, np.vstack((f, f+len(v))))
    source = np.array([0., 0., 400.])
    rays = vertices[faces].mean(axis=1)-source
    front = np.flatnonzero(np.einsum("ij,ij->i", normals, -rays) > 1e-9)
    visible, _ = visible_faces(vertices, faces, normals, source, front)
    assert np.all(visible[front < len(f)])          # the near sphere is seen
    assert not np.any(visible[front >= len(f)])     # the far one is fully hidden behind it
    # A source-facing triangle that reaches the focal plane disables the kernel.
    assert visible_faces(vertices, faces, normals, np.array([0., 0., 105.]), front[:3]) is None
