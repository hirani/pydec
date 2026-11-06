"""The vectorised whitney_innerproduct against the per-simplex reference loop and
against the per-simplex barycentric_gradients.
"""
import os
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.spatial import Delaunay

from pydec.dec import simplicial_complex
from pydec.fem.innerproduct import (whitney_innerproduct, _whitney_innerproduct_reference,
                                    barycentric_gradients, _barycentric_gradients_batch)
from pydec.io import read_mesh

mesh_dir = os.path.join(os.path.split(__file__)[0], '..', '..', 'dec', 'tests', 'meshes')


def example_meshes():
    rng = np.random.default_rng(99)
    M = {}
    x = np.sort(rng.random(20))
    T = np.array([[i, i + 1] for i in range(19)]); T[::3] = T[::3, ::-1]
    M['segments_1d'] = (x[:, None], T)
    th = np.linspace(0, 3, 20)
    M['polyline_2d'] = (np.c_[np.cos(th), np.sin(th)], T)
    V = rng.random((120, 2)); T = Delaunay(V).simplices
    M['tri_delaunay'] = (V, T)
    V = rng.random((150, 3)); T = Delaunay(V).simplices
    M['tet_delaunay'] = (V, T)
    Ts = T.copy()
    for row in Ts:
        rng.shuffle(row)
    M['tet_shuffled'] = (V, Ts)
    V = rng.random((40, 4)); T = Delaunay(V).simplices
    M['tet_4d'] = (V, T)
    for name in ['sphere3', 'cube', 'torus_tet']:
        m = read_mesh(os.path.join(mesh_dir, name + '.xml'))
        M[name] = (np.asarray(m['vertices'], dtype=float), np.asarray(m['elements']))
    return M


def test_gradients_batch_matches_single():
    rng = np.random.default_rng(5)
    for K in [1, 2, 3, 4]:
        for M in [1, 2, K + 1]:
            if M > K + 1: continue
            pts = rng.random((30, M, K))
            g = _barycentric_gradients_batch(pts)
            for i in range(len(pts)):
                assert_allclose(g[i], barycentric_gradients(pts[i]), rtol=1e-11, atol=1e-11)


def test_whitney_matches_reference():
    """same sparsity structure; entries within eps times the conditioning of the
    worst element, which for random Delaunay meshes with slivers is about 1e-9"""
    for name, (V, T) in example_meshes().items():
        sc = simplicial_complex((V, T))
        tol = 1e-12 if name in ('segments_1d', 'polyline_2d', 'sphere3', 'cube', 'torus_tet') else 1e-7
        for k in range(sc.complex_dimension() + 1):
            A = whitney_innerproduct(sc, k).tocsr(); A.sort_indices()
            B = _whitney_innerproduct_reference(sc, k).tocsr(); B.sort_indices()
            assert_array_equal(A.indptr, B.indptr, err_msg=f'{name} k={k}')
            assert_array_equal(A.indices, B.indices, err_msg=f'{name} k={k}')
            scale = np.abs(B.data).max()
            assert_allclose(A.data, B.data, rtol=0, atol=tol * scale, err_msg=f'{name} k={k}')


def test_whitney_symmetric_and_positive_on_oriented_mesh():
    m = read_mesh(os.path.join(mesh_dir, 'torus_tet.xml'))
    V, T = np.asarray(m['vertices'], dtype=float), np.asarray(m['elements'])
    sc = simplicial_complex((V, T))
    assert (sc[3].primal_volume > 0).all()
    for k in range(4):
        A = whitney_innerproduct(sc, k)
        assert abs(A - A.T).max() < 1e-14 * abs(A).max()
        assert np.linalg.eigvalsh(A.toarray()).min() > 0
