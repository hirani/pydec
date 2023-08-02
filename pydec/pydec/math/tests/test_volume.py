from pydec.testing import *

from scipy import fabs, random, rand, array, sqrt

from pydec.math.volume import unsigned_volume, signed_volume


def test_unsigned_volume():
    cases = []
    cases.append((array([[1]]), 1))
    cases.append((array([[1],[10]]), 9))
    cases.append((array([[0,0],[1,1]]),  sqrt(2)))
    cases.append((array([[0,0],[0,1],[1,0]]), 1.0/2.0))
    cases.append((array([[0,0],[0,1],[1,0]]), 1.0/2.0))
    cases.append((array([[5,5],[5,6],[6,5]]), 1.0/2.0))
    cases.append((array([[0,0,0],[0,0,1],[0,1,0]]), 1.0/2.0))
    
    for s,v in cases:
        assert_almost_equal(unsigned_volume(s), v)

def test_signed_volume():
    cases = []
    cases.append((array([[1],[2]]), 1))
    cases.append((array([[5.5],[-10]]), -15.5))
    cases.append((array([[0,0],[1,1],[1,0]]), -1.0/2.0))
    cases.append((array([[0,0],[1,0],[1,1]]),  1.0/2.0))
    cases.append((array([[0,0],[0,1],[1,0]]), -1.0/2.0))        
    cases.append((array([[5,5],[5,6],[6,5]]), -1.0/2.0))
    cases.append((array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]), 1.0/6.0))
    
    for s,v in cases:
        assert_almost_equal(signed_volume(s), v)
        
def test_both():
    """signed and unsigned volumes should agree up to sign"""

    random.seed(0) #make tests repeatable                 
    for N in range(1,10):
        pts = rand(N+1,N)
        assert_almost_equal(fabs(signed_volume(pts)), unsigned_volume(pts))
      
