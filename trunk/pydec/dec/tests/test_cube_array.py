from pydec.testing import *

from scipy import random, arange, alltrue, array

from pydec.dec.cube_array import cube_array_search

class test_cube_array_search(TestCase):
    def setUp(self):
        random.seed(0)


    def test_simple1(self):
        k_face_array = array([[0,0,0],
                              [0,0,1],
                              [0,1,0],
                              [1,0,1]])

        k_faces = array([[0,1,0],
                         [0,0,1]])

        indices = cube_array_search(k_face_array, k_faces)

        assert_equal(indices, array([2,1]))

    def test_simple2(self):
        k_face_array = array([[0,1,0],
                              [0,0,2],
                              [0,2,0],
                              [1,0,3]])

        k_faces = array([[0,1,0],
                         [0,2,0],
                         [0,0,2],
                         [0,0,2],
                         [1,0,3],
                         [0,1,0]])

        indices = cube_array_search(k_face_array, k_faces)

        assert_equal(indices, array([0,2,1,1,3,0]))

