"""The vectorised geometry of simplicial_complex (circumcenters, primal and
dual volumes) against the per-simplex functions in pydec.math and the
per-simplex dual volume recursion kept as simplicial_complex._compute_dual_volume_reference.
"""
import os
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.spatial import Delaunay

from pydec.dec import simplicial_complex
from pydec.math import (circumcenter, circumcenter_barycentric, circumcenters,
                        unsigned_volume, signed_volume, unsigned_volumes, signed_volumes)
from pydec.io import read_mesh

mesh_dir = os.path.join(os.path.split(__file__)[0], 'meshes')


def perturbed(V, T, rng, frac):
    """move the vertices by frac of the mean edge length: no longer Delaunay, obtuse elements"""
    h = np.linalg.norm(V[T[:, 1]] - V[T[:, 0]], axis=1).mean()
    return V + rng.normal(size=V.shape) * frac * h


def shuffled_columns(T, rng):
    """the same simplices with their vertices listed in random order"""
    T = T.copy()
    for row in T:
        rng.shuffle(row)
    return T


def grid_tets(n):
    """Kuhn subdivision of an n x n x n cube grid: cospherical vertices, so
    some barycentric circumcenter coordinates are exactly zero"""
    ax = np.linspace(0, 1, n + 1)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij')
    V = np.c_[X.ravel(), Y.ravel(), Z.ravel()]
    idx = lambda i, j, k: (i * (n + 1) + j) * (n + 1) + k
    T = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for p in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]:
                    c = np.array([i, j, k]); tet = [idx(*c)]
                    for a in p:
                        c = c.copy(); c[a] += 1; tet.append(idx(*c))
                    T.append(tet)
    return V, np.array(T)


def example_meshes():
    rng = np.random.default_rng(2026)
    M = {}
    x = np.sort(rng.random(20))
    T = np.array([[i, i + 1] for i in range(19)]); T[::3] = T[::3, ::-1]
    M['segments_1d'] = (x[:, None], T)
    th = np.linspace(0, 3, 20)
    M['polyline_2d'] = (np.c_[np.cos(th), np.sin(th)], T)
    V = rng.random((120, 2)); T = Delaunay(V).simplices
    M['tri_delaunay'] = (V, T)
    M['tri_perturbed'] = (perturbed(V, T, rng, 0.3), T)
    M['tri_shuffled'] = (V, shuffled_columns(T, rng))
    V = rng.random((150, 3)); T = Delaunay(V).simplices
    M['tet_delaunay'] = (V, T)
    M['tet_perturbed'] = (perturbed(V, T, rng, 0.3), T)
    M['tet_shuffled'] = (V, shuffled_columns(T, rng))
    M['tet_grid'] = grid_tets(3)
    V = rng.random((40, 4)); T = Delaunay(V).simplices
    M['tet_4d'] = (V, T)
    for name in ['sphere3', 'cube', 'torus_tet']:
        m = read_mesh(os.path.join(mesh_dir, name + '.xml'))
        M[name] = (np.asarray(m['vertices'], dtype=float), np.asarray(m['elements']))
    return M


def test_batch_functions_match_single():
    rng = np.random.default_rng(7)
    for K in [1, 2, 3, 4]:
        for M in range(1, K + 2):
            pts = rng.random((25, M, K))
            c, bc = circumcenters(pts)
            for i in range(len(pts)):
                assert_allclose(bc[i], circumcenter_barycentric(pts[i]), rtol=1e-12, atol=1e-12)
                assert_allclose(c[i], circumcenter(pts[i])[0], rtol=1e-12, atol=1e-12)
                assert_allclose(unsigned_volumes(pts)[i], unsigned_volume(pts[i]), rtol=1e-12)
            if M == K + 1:
                for i in range(len(pts)):
                    assert_allclose(signed_volumes(pts)[i], signed_volume(pts[i]), rtol=1e-12)


def test_geometry_matches_per_simplex():
    for name, (V, T) in example_meshes().items():
        sc = simplicial_complex((V, T))
        n, N = sc.complex_dimension(), sc.embedding_dimension()
        for k, data in enumerate(sc):
            for i, s in enumerate(data.simplices):
                pts = V[s, :]
                assert_allclose(data.circumcenter[i], circumcenter(pts)[0], rtol=1e-12, atol=1e-12, err_msg=name)
                assert_allclose(data.bary_circumcenter[i], circumcenter_barycentric(pts), rtol=1e-12, atol=1e-12, err_msg=name)
            if k == N:
                expected = [signed_volume(V[s, :]) for s in sc.simplices]
            else:
                expected = [unsigned_volume(V[s, :]) for s in data.simplices]
            assert_allclose(data.primal_volume, expected, rtol=1e-12, atol=1e-15, err_msg=name)


def test_dual_volume_matches_reference():
    for name, (V, T) in example_meshes().items():
        sc = simplicial_complex((V, T))
        fast = [data.dual_volume.copy() for data in sc]
        sc._compute_dual_volume_reference()
        for k, data in enumerate(sc):
            scale = np.abs(data.dual_volume).max()
            assert_allclose(fast[k], data.dual_volume, rtol=0, atol=1e-13 * max(scale, 1.0), err_msg=f'{name} dim {k}')
        assert_array_equal(fast[-1], 1.0)


def test_dual_volume_identity():
    """sum over k-simplices of |s| |*s| = C(n,k) |K| for a flat n-complex"""
    from math import comb
    for name, (V, T) in example_meshes().items():
        sc = simplicial_complex((V, T))
        n, N = sc.complex_dimension(), sc.embedding_dimension()
        if n != N:
            continue
        total = np.abs(sc[n].primal_volume).sum()
        for k, data in enumerate(sc):
            # 1e-8, not 1e-12: the circumcenter linear systems are ill conditioned
            # for the cospherical grid tets, and the identity sums many signed terms
            assert_allclose(np.dot(np.abs(data.primal_volume), data.dual_volume), comb(n, k) * total,
                            rtol=1e-8, err_msg=f'{name} dim {k}')
