import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from skindose_verifier.surface_views import build_views, save_surface_views


def _cylinder(radius=100., length=300., n_theta=72, n_z=30):
    """Closed, outward-oriented cylinder along +z (CT LPS: z craneal, y posterior)."""
    theta = np.linspace(0, 2*np.pi, n_theta, endpoint=False)
    z = np.linspace(0, length, n_z)
    ring = np.column_stack([radius*np.cos(theta), radius*np.sin(theta)])
    vertices = np.array([[x, y, zz] for zz in z for x, y in ring])
    faces = []
    for k in range(n_z - 1):
        for i in range(n_theta):
            a, b = k*n_theta + i, k*n_theta + (i + 1) % n_theta
            c, d = a + n_theta, b + n_theta
            faces += [[a, b, d], [a, d, c]]
    bottom, top = len(vertices), len(vertices) + 1
    vertices = np.vstack([vertices, [0, 0, 0], [0, 0, length]])
    for i in range(n_theta):
        faces.append([bottom, (i + 1) % n_theta, i])
        faces.append([top, (n_z - 1)*n_theta + i, (n_z - 1)*n_theta + (i + 1) % n_theta])
    faces = np.asarray(faces)
    tri = vertices[faces]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    return vertices, faces, normals


def test_posterior_field_is_seen_from_behind_and_hidden_from_the_front(tmp_path):
    vertices, faces, normals = _cylinder()
    centroids = vertices[faces].mean(axis=1)
    assert np.all(np.einsum("ij,ij->i", normals, centroids - [0, 0, 150]) > 0)  # outward
    posterior_field = (normals[:, 1] > .7) & (np.abs(centroids[:, 2] - 150) < 40)
    dose = np.where(posterior_field, .4, 0.)
    views = {view.name: view for view in build_views(vertices, faces, ct_lps=True)}
    seen_behind = np.unique(views["Posterior"].image_face[views["Posterior"].image_face >= 0])
    seen_front = np.unique(views["Anterior"].image_face[views["Anterior"].image_face >= 0])
    assert np.isin(np.flatnonzero(posterior_field), seen_behind).mean() > .9
    assert not np.isin(np.flatnonzero(posterior_field), seen_front).any()
    # Seen from behind, the patient's left (+x in LPS) is on the viewer's left.
    assert (views["Posterior"].left_label, views["Posterior"].right_label) == ("I", "D")
    assert (views["Anterior"].left_label, views["Anterior"].right_label) == ("D", "I")
    assert (views["Lateral izquierda"].left_label, views["Lateral izquierda"].right_label) == ("A", "P")
    # No pinholes inside the silhouette.
    image = views["Posterior"].image_face
    rows = image[image.shape[0]//4: 3*image.shape[0]//4, image.shape[1]//3: 2*image.shape[1]//3]
    assert np.all(rows >= 0)
    palette = LinearSegmentedColormap.from_list("d", ["#2ca25f", "#ffe34d", "#fb8c00", "#d7191c"])
    path = tmp_path/"vistas.png"
    save_surface_views(path, vertices, faces, dose, normals, ct_lps=True, palette=palette)
    assert path.stat().st_size > 10_000
