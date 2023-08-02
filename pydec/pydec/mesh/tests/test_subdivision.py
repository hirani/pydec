#import unittest
#
#from pydec.mesh import cube_grid,triangulate_ncube
#from pydec.math.volume import unsigned_volume,signed_volume
#from numpy import allclose
#
#
#class TestSudivision(unittest.TestCase):
#    def setUp(self):
#        pass
#
#    def test_triangulate_ncube(self):
#        for n in range(1,6):
#            v,i = cube_grid((1,)*n)
#            v,i = triangulate_ncube(v,i)
#            
#            unsigned_sum = 0
#            signed_sum = 0
#            for row in i:
#                unsigned_sum += unsigned_volume(v[row])
#                signed_sum   += signed_volume(v[row])
#                
#            self.assert_(allclose(unsigned_sum,1))
#            self.assert_(allclose(signed_sum,1))
#
#
