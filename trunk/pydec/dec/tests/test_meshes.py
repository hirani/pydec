from pydec.testing import *

import os

from scipy import random, arange, alltrue, array, rand, allclose, \
        concatenate, zeros, shape, sparse, matrix, sqrt

from pydec.dec import simplicial_complex
from pydec.mesh import simplicial_mesh,simplex
from pydec.io import read_mesh



base_dir = os.path.split(__file__)[0]
mesh_dir = os.path.join( base_dir, 'meshes')

mesh_filenames = ['cube.xml','sphere3.xml','torus_tet.xml']

meshes = {}
for name in mesh_filenames:
    meshes[name] = read_mesh(os.path.join(mesh_dir,name))


         
def test_exactness():
    for mesh in meshes.itervalues():
        sc = simplicial_complex(mesh)

        cmplx = sc.chain_complex()

        for B1,B2 in zip(cmplx[:-1],cmplx[1:]):
            assert_equal((B1*B2).nnz,0)

def test_faces():
    """check whether faces are correct"""
    for mesh in meshes.itervalues():
        sc = simplicial_complex(mesh)

        for face_data,cell_data in zip(sc[:-1],sc[1:]):
            #compute all faces and store in a set

            boundary_set = set()
            for row in cell_data.simplices:
                boundary_set.update(simplex(row).boundary())

            face_set = set([simplex(row) for row in face_data.simplices])

            assert_equal(boundary_set, face_set)

            
def test_boundary():
    """check boundary operators (assumes correct faces)"""
    
    for mesh in meshes.itervalues():
        sc = simplicial_complex(mesh)
        assert_equal(sc[0].boundary.shape,(1,sc[0].simplices.shape[0]))
        assert_equal(sc[0].boundary.nnz,0)

        for face_data,cell_data in zip(sc[:-1],sc[1:]):
            B = cell_data.boundary
           
            assert_equal(B.shape,(face_data.num_simplices,cell_data.num_simplices))
            assert_equal(B.nnz,cell_data.num_simplices*(cell_data.dim+1))
            
            for row,parity in zip(cell_data.simplices,cell_data.simplex_parity):
                s = simplex(row)
                
                s_index = cell_data.simplex_to_index[s]
                
                for f in s.boundary():
                    f_index = face_data.simplex_to_index[f]

                    assert_equal(B[f_index,s_index],(-1)**(f.parity ^ int(parity)))
                
