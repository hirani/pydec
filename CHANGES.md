# Changes

## Unreleased (branch `fast-hodge`)

### Circumcentric Hodge star assembly: minutes to seconds

`simplicial_complex.compute_dual_volume`, `compute_circumcenters`,
`compute_bary_circumcenters` and `compute_primal_volume` now work on all simplices of a
dimension at once. Before, the dual volumes were computed by a Python recursion over every
chain of faces of every top simplex (41 calls per tetrahedron, each with a dictionary
lookup and a `scipy.linalg.det` on a 1x1 or 2x2 matrix), which took 52–71 minutes for the
two Hodge stars of a 1.07M-tetrahedron mesh. The same mesh now takes 20 s; a 13k-tet mesh
0.2 s instead of 9.4 s.

Nothing visible changes: the attributes `circumcenter`, `bary_circumcenter`,
`primal_volume`, `dual_volume`, `star` and `star_inv` have the same shapes, the same
indexing and the same sign convention (the sign of a dual simplex is the product, along
its chain of faces, of the sign of the barycentric coordinate of each circumcenter at the
vertex opposite the face; Hirani, Kalyanaraman and VanderZee, *Delaunay Hodge star*,
Computer-Aided Design 45 (2013)). The new code was checked against the old on 1-D to 4-D
meshes, Delaunay and not, with shuffled vertex order, with cospherical vertices, and on
surfaces in 3-D, and on several hundred meshes from the authors' projects: circumcenters
are bitwise identical and dual volumes agree to 1e-13 of their largest value.

The old recursion is kept as `simplicial_complex._compute_dual_volume_reference` for
tests; `dec/tests/test_dual_volume.py` compares the two.

One difference, on meshes that list the same top simplex twice (an error in the mesh): the
old code gave one copy a dual volume of 2 and the other 0, an accident of its dictionary
lookup; the new code gives each copy 1. The dual volumes of the faces are the same in
both, and both count the faces of the duplicated simplex twice.

New functions in `pydec.math`, each acting on a stack of simplices: `circumcenters`,
`circumcenters_barycentric`, `unsigned_volumes`, `signed_volumes`. The per-simplex
`circumcenter`, `circumcenter_barycentric`, `unsigned_volume` and `signed_volume` are
unchanged apart from `numpy.bmat` becoming `numpy.block`, which removes a
`PendingDeprecationWarning` and does not change the numbers.

### `star_inv` is built on its own and warns about zero dual volumes

`construct_hodge` used to build `star` and `star_inv` together, so asking for `star` on a
mesh with a zero dual volume produced NumPy's divide-by-zero warning. Now `sc[k].star`
calls `construct_hodge_star` and `sc[k].star_inv` calls `construct_hodge_star_inv`;
`construct_hodge` still builds both. The values and the sign convention,
`star_inv * star = (-1)^(k(n-k))`, are unchanged. When `star_inv` is built, a
`RuntimeWarning` says how many dual volumes of each dimension are zero (infinite entries)
or below `tol` (default 1e-12) of the largest (huge entries); nothing is clipped or
regularised. Zero dual volumes occur on meshes with cospherical vertices, for instance a
cube cut into six tetrahedra that share the cube's centre, and tiny ones by cancellation
on meshes that are not Delaunay.

### NumPy 2.4 compatibility

`whitney_innerproduct(sc, 1)` raised `ValueError: setting an array element with a
sequence` on NumPy >= 2.4, because for 1-forms it stored a 1x1 array into a float
(deprecated in NumPy 1.25). Fixed in `fem/innerproduct.py`.

### Test suite

The tests did not collect under pytest with current SciPy, which no longer re-exports
NumPy names, and `pydec.testing` exported only `TestCase`. The tests now collect and run:
`python -m pytest pydec` from the repository top with the package on `PYTHONPATH`, or
after `pip install`. Python 3 and NumPy 2 bugs in `cube_array.py`, `io/arrayio.py` and
`math/parity.py` that the tests exposed were fixed; no assertion was changed.
