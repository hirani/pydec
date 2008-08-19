from pydec.testing import *

from scipy import arange, array

from pydec.dec.cube_array import * 

class test_cube_array_search(TestCase):
    def test_simple1(self):
        k_face_array = array([[0,0,0],
                              [0,0,1],
                              [0,1,0],
                              [1,0,1]])
        rows = array([0,2,1,1,3,0])
        k_faces = k_face_array[rows, :]
        assert_equal(cube_array_search(k_face_array, k_faces), rows)

    def test_simple2(self):
        k_face_array = array([[0,1,0],
                              [0,0,2],
                              [0,2,0],
                              [1,0,3]])
        rows = array([0,2,1,1,3,0])
        k_faces = k_face_array[rows, :]
        assert_equal(cube_array_search(k_face_array, k_faces), rows)
    
    def test_simple3(self):
        k_face_array = array([[0,1,0,1],
                              [0,0,1,2],
                              [0,2,0,2],
                              [1,0,0,1]])

        rows = array([2,1,0,0,3,1,2])
        k_faces = k_face_array[rows,:]
        assert_equal(cube_array_search(k_face_array, k_faces), rows)


class test_cube_array_boundary(TestCase):
    def test_boundary1_1D(self):
        """boundary_1 of 1D mesh [XXX]"""
        cubes = array([[0,0],
                       [1,0],
                       [2,0]])
        dim = 1

        faces,boundary = cube_array_boundary(cubes, dim)

        assert_equal(faces, [[0],
                             [1],
                             [2],
                             [3]] )
        assert_equal(boundary.toarray(), [[-1, 0, 0],
                                          [ 1,-1, 0],
                                          [ 0, 1,-1],
                                          [ 0, 0, 1]] )
    def test_boundary1_1D_with_holes(self):
        """boundary_1 of 1D mesh [X00X0X]"""
        cubes = array([[0,0],
                       [3,0],
                       [5,0]])
        dim = 1

        faces,boundary = cube_array_boundary(cubes, dim)

        assert_equal(faces, [[0],
                             [1],
                             [3],
                             [4],
                             [5],
                             [6]] )
        assert_equal(boundary.toarray(), [[-1, 0, 0],
                                          [ 1, 0, 0],
                                          [ 0,-1, 0],
                                          [ 0, 1, 0],
                                          [ 0, 0,-1],
                                          [ 0, 0, 1]] )
    
    def test_boundary1_2D(self):
        """boundary_1 of 2D mesh [[X]]"""
        cubes = array([[0,0,0],
                       [0,0,1],
                       [0,1,0],
                       [1,0,1]])
        dim = 1

        faces,boundary = cube_array_boundary(cubes, dim)

        assert_equal(faces, [[0,0],
                             [0,1],
                             [1,0],
                             [1,1]] )
        assert_equal(boundary.toarray(), [[-1,-1, 0, 0],
                                          [ 0, 1,-1, 0],
                                          [ 1, 0, 0,-1],
                                          [ 0, 0, 1, 1]] )


    def test_boundary2_2D(self):
        """boundary_2 of 2D mesh [[XX],[XX]]"""
        cubes = array([[0,0,0,1],
                       [0,1,0,1],
                       [1,0,0,1],
                       [1,1,0,1]])
        dim = 2

        faces,boundary = cube_array_boundary(cubes, dim)

        assert_equal(faces, [[0,0,0],
                             [0,0,1],
                             [0,1,0],
                             [0,1,1],
                             [0,2,0],
                             [1,0,0],  
                             [1,0,1],
                             [1,1,0],
                             [1,1,1],
                             [1,2,0],
                             [2,0,1],
                             [2,1,1]] )
        assert_equal(boundary.toarray(), [[ 1, 0, 0, 0],
                                          [-1, 0, 0, 0],
                                          [-1, 1, 0, 0],
                                          [ 0,-1, 0, 0],
                                          [ 0,-1, 0, 0],
                                          [ 0, 0, 1, 0],
                                          [ 1, 0,-1, 0],
                                          [ 0, 0,-1, 1],
                                          [ 0, 1, 0,-1],
                                          [ 0, 0, 0,-1],
                                          [ 0, 0, 1, 0],
                                          [ 0, 0, 0, 1]] )

    def test_boundary2_2D_with_holes(self):
        """boundary_2 of 2D mesh [[X0X]]"""
        cubes = array([[0,0,0,1],
                       [2,0,0,1]])
        dim = 2

        faces,boundary = cube_array_boundary(cubes, dim)

        assert_equal(faces, [[0,0,0],
                             [0,0,1],
                             [0,1,0],
                             [1,0,1],
                             [2,0,0], 
                             [2,0,1],
                             [2,1,0],
                             [3,0,1]])
        assert_equal(boundary.toarray(), [[ 1, 0],
                                          [-1, 0],
                                          [-1, 0],
                                          [ 1, 0],
                                          [ 0, 1],
                                          [ 0,-1],
                                          [ 0,-1],
                                          [ 0, 1]] )



    def test_boundary2_3D(self):
        """boundary_2 of a 3D cube"""
        cubes = array([[0,0,0,0,1,2]])
        dim = 3

        faces,boundary = cube_array_boundary(cubes, dim)

        assert_equal(faces, [[0,0,0,0,1],
                             [0,0,0,0,2],
                             [0,0,0,1,2],
                             [0,0,1,0,1],
                             [0,1,0,0,2],
                             [1,0,0,1,2]] )
        assert_equal(boundary.toarray(), [[-1],
                                          [ 1],
                                          [-1],
                                          [ 1],
                                          [-1],
                                          [ 1]] )


