from pydec.testing import *

from scipy import random,rand,array,reshape,sqrt,sum,allclose

from pydec.math.circumcenter import circumcenter,is_wellcentered

class TestCircumcenter(TestCase):
    def setUp(self):	
        random.seed(0) #make tests repeatable                 

    def test_is_well_centered(self):
        self.assert_(is_wellcentered(array([[1]])))
        self.assert_(is_wellcentered(array([[0,1]])))
        self.assert_(is_wellcentered(array([[0],[1]])))
        self.assert_(is_wellcentered(array([[0,1],[1,0]])))
        self.assert_(is_wellcentered(array([[0,0],[1,0],[0.5,sqrt(3)/2]])))
        
    
    def test_spheres(self):
        """points on an N dimensional sphere with known center and radius"""
        for N in range(2,10):            
            pts = rand(N+1,N)            
            true_center = rand(N)
            true_radius = 1+10*rand()
            pts = pts/reshape(sqrt(sum(pts*pts,axis=1)),(N+1,1))
            pts *= true_radius
            pts += true_center
            (center,radius) = circumcenter(pts)
            assert_almost_equal(center,true_center)
            assert_almost_equal(radius,true_radius)
            
    def test_random(self):
        """M random points in N dimensional space (where M <= N + 1)"""        
        for N in range(1,10):
            for M in range(1,N+1):
                pts = rand(M,N)
                (center,radius) = circumcenter(pts)           
                distances = sqrt(sum((pts - center)**2,axis=1))
                assert_almost_equal(max(distances),min(distances))
                             

