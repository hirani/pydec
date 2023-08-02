from pydec.testing import *

from scipy import array

from pydec.mesh.simplex import simplex,simplicial_mesh

  
class test_simplicial_mesh(TestCase):
    def setUp(self):
        self.meshes = []
        self.meshes.append((array([[0],[1]]),array([[0,1]])))
        self.meshes.append((array([[0],[1],[2]]),array([[0,1],[1,2]])))
        self.meshes.append((array([[0],[1],[2],[3]]),array([[0,1],[1,2],[2,3]])))
        self.meshes.append((array([[0,0],[1,0],[0,1]]),array([[0,1],[1,2],[2,0]])))        
        self.meshes.append((array([[0,0],[1,0],[0,1]]),array([[0,1,2]])))        
        self.meshes.append((array([[0,0],[1,0],[0,1],[1,1]]),array([[0,1,2],[1,3,2]])))
        self.meshes.append((array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]),array([[0,1,2,3]])))
        self.meshes.append((array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1]]),array([[0,1,2,3],[1,2,3,4]])))
                

    def check_boundary(self):
        """Test the boundary() method"""
        
        boundaries = []
        boundaries.append([simplex([0],parity=1),simplex([1],parity=0)])
        boundaries.append([simplex([0],parity=1),simplex([2],parity=0)])
        boundaries.append([simplex([0],parity=1),simplex([3],parity=0)])
        boundaries.append([])
        boundaries.append([simplex([0,1]),simplex([1,2]),simplex([2,0])])
        boundaries.append([simplex([0,1]),simplex([1,3]),simplex([3,2]),simplex([2,0])])
        boundaries.append([simplex([1,0,2]),simplex([0,1,3]),simplex([2,0,3]),simplex([1,2,3])])
        boundaries.append([simplex([1,0,2]),simplex([0,1,3]),simplex([2,0,3]),simplex([2,3,4]),simplex([1,4,3]),simplex([1,2,4])])
        
        cases = zip(self.meshes,boundaries)
        
        for (v,s),b in cases:
            sm = simplicial_mesh(v,s)                           
            sb = sm.boundary()
            assert_equal(sb,set(b))
            
            for x in sb:         
                assert_equal(x.parity,b[b.index(x)].parity)
                
                
            
    def check_skeleton(self):
        """Test the skeleton() method"""

        assert_equal(simplicial_mesh(array([[0]]),array([[0]])).skeleton(0), \
                        set([simplex([0])]))        
        assert_equal(simplicial_mesh(array([[0],[1]]),array([[0,1]])).skeleton(0), \
                        set([simplex([0]),simplex([1])]))        
        assert_equal(simplicial_mesh(array([[0],[1]]),array([[0,1]])).skeleton(1), \
                        set([simplex([0,1])]))
        assert_equal(simplicial_mesh(array([[0],[1],[2]]),array([[0,1],[1,2]])).skeleton(0),\
                        set([simplex([0]),simplex([1]),simplex([2])]))
        assert_equal(simplicial_mesh(array([[0],[1],[2]]),array([[0,1],[1,2]])).skeleton(1),\
                        set([simplex([0,1]),simplex([1,2])]))
        assert_equal(simplicial_mesh(array([[0,0],[1,0],[0,1]]),array([[0,1,2]])).skeleton(1),\
                        set([simplex([0,1]),simplex([1,2]),simplex([2,0])]))
        
        
