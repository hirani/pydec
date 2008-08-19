from pydec.testing import *

from scipy import array, alltrue, array, sqrt, sparse

from pydec.dec import simplicial_complex
from pydec.dec.cochain import d, laplace_beltrami, laplace_derham, star
from pydec.mesh.simplex import simplex


all_cases = []
#unit line segments
all_cases.append((array([[0],[1]]),array([[0,1]])))
all_cases.append((array([[1],[0]]),array([[1,0]])))        
#regular triangle, edge length 1
all_cases.append((array([[0,0],[1,0],[0.5,sqrt(3.0)/2.0]]),array([[0,1,2]])))        
#regular tet, edge length sqrt(2)
all_cases.append((array([[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 0]]),array([[0,1,2,3]])))               
# 4 triangles in a square
all_cases.append((array([[0,0],[1,0],[1,1],[0,1],[0.5,0.5]]),array([[0,1,4],[1,2,4],[2,3,4],[0,4,3]])))
#two line segments
all_cases.append((array([[0],[1],[2]]),array([[0,1],[1,2]])))
#three line segments
all_cases.append((array([[0],[1],[2],[3]]),array([[0,1],[1,2],[2,3]])))



def test_three_edges():
    V = array([[0],[1],[2],[3]])  
    S = array([[0,1],[1,2],[2,3]]) 
    sc = simplicial_complex((V,S))
    f0 = sc.get_cochain(0)
    
    assert_equal(f0.complex,sc)
    assert_equal(f0.k,0)
    assert_equal(f0.is_primal,True)
    assert_equal(f0.v,array([0,0,0,0]))
            
    f0.v[:] = 10
    f1 = d(f0)
   
    assert_equal(f1.complex,sc)
    assert_equal(f1.k,1)
    assert_equal(f1.is_primal,True)
    assert_equal(f1.v,array([0,0,0]))

def test_get_set():
    V = array([[0,0],[1,0],[0.5,sqrt(3.0)/2.0]]) 
    S = array([[0,1,2]]) 
    sc = simplicial_complex((V,S))
    
    c0 = sc.get_cochain(0)
    assert_equal(c0[0],0)
    assert_equal(c0[0],0)
    assert_equal(c0[0],0)        
    
    c0[simplex([0],parity=0)] = 10
    c0[simplex([1],parity=0)] = 20
    c0[simplex([2],parity=1)] = 30        
    assert_equal(c0[simplex([0],parity=0)],10)
    assert_equal(c0[simplex([1],parity=1)],-20)
    assert_equal(c0[simplex([2],parity=1)],30)
    assert_equal(c0[0],10)
    assert_equal(c0[1],20)
    assert_equal(c0[2],-30)
        
        
class TestCochainFunctions(TestCase):
    def test_d(self):
        for case in all_cases:
            sc = simplicial_complex(case)
            N  = sc.complex_dimension()
            for p in range(sc.complex_dimension()):
                df = d(sc.get_cochain_basis(p))
                assert_equal( (df.v - sc[p].d).nnz, 0)
                
                # test exactness
                assert_equal(d(df).v.nnz,0)

    def test_dual_d(self):
        for case in all_cases:
            sc = simplicial_complex(case)
            N  = sc.complex_dimension()
            for p in range(sc.complex_dimension()):
                df = d(sc.get_cochain_basis(p, is_primal=False))
                assert_equal( (df.v - (-1)**p*sc[N-p].boundary).nnz, 0)
                
                # test exactness
                assert_equal(d(df).v.nnz,0)

        
    def test_star(self):
        for case in all_cases:
            sc = simplicial_complex(case)
            N  = sc.complex_dimension()
            for p in range(N):
                if not (sc[p].dual_volume > 0).all(): continue #skip non-wellcentered

                result   = star(star(sc.get_cochain_basis(p))).v
                expected = (-1)**(p*(N - p)) * sparse.identity(sc[p].num_simplices)
                assert_almost_equal(result.todense(),expected.todense())
        
        
    def test_delta(self):    
        pass
        
    def test_laplacian(self):
        """Check that Laplace-Beltrami and Laplace-DeRham agree for 0-forms
        
        This makes sure that d() and star() work in the -1 and N+1 cases
        """
        for case in all_cases:
            sc = simplicial_complex(case)
            f  = sc.get_cochain_basis(0)
            assert_equal((laplace_beltrami(f) - laplace_derham(f)).v.nnz, 0)        
        

