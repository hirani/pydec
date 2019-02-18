from pydec.testing import *

from scipy import array,array,sparse,zeros,ones,eye,allclose,alltrue, \
                isreal,real,dot,concatenate,sqrt,shape,ix_
from scipy.special import factorial, comb
from scipy.linalg import eigvals,inv,det

from pydec.dec import simplicial_complex, regular_cube_complex
from pydec.mesh import simplex, regular_cube_mesh
from pydec.math.combinatorial import combinations
from pydec.fem.innerproduct import whitney_innerproduct, barycentric_gradients, massmatrix_rowcols, \
     regular_cube_innerproduct
    

class TestGradients(TestCase):    
    def test_gradients(self):
        cases = []
        cases.append((array([[0],[1]]),array([[-1],[1]])))
        cases.append((array([[0,1,2],[1,2,3]]),array([[-1.0/3.0,-1.0/3.0,-1.0/3.0],[1.0/3.0,1.0/3.0,1.0/3.0]])))
        cases.append((array([[0,0],[1,0],[0,1]]),array([[-1,-1],[1,0],[0,1]])))
        cases.append((array([[0,0],[1,0],[0.5,sqrt(3.0)/2.0]]), \
                    array([[-1.0,-1.0/sqrt(3.0)],[1.0,-1.0/sqrt(3.0)],[0.0,2.0/sqrt(3.0)]])))
        cases.append((array([[0,0,0],[1,0,0],[0,1,0]]),array([[-1,-1,0],[1,0,0],[0,1,0]])))        
        cases.append((array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]), array([[-1,-1,-1],[1,0,0],[0,1,0],[0,0,1]])))        
        cases.append((array([[0,0,0,0],[1,0,0,0],[0,1,0,0],[0,0,1,0]]), \
                    array([[-1,-1,-1,0],[1,0,0,0],[0,1,0,0],[0,0,1,0]])))        
    
        for pts,grads in cases:                        
            assert_almost_equal(grads,barycentric_gradients(pts))

class TestMassmatrixRowCols(TestCase):    
    def test_simple(self):
        V = array([[0,0],[1,0],[0,1]])
        S = array([[0,1,2]])

        sc = simplicial_complex((V,S))

        rows,cols = massmatrix_rowcols(sc,0)

        assert_equal(rows,array([0,0,0,1,1,1,2,2,2]))
        assert_equal(cols,array([0,1,2,0,1,2,0,1,2]))


        rows,cols = massmatrix_rowcols(sc,1)

        edges = [simplex(x) for x in combinations(range(3),2)]
        edge0 = sc[1].simplex_to_index[edges[0]]
        edge1 = sc[1].simplex_to_index[edges[1]]
        edge2 = sc[1].simplex_to_index[edges[2]]

        assert_equal(rows,array([edge0,edge0,edge0,edge1,edge1,edge1,edge2,edge2,edge2]))
        assert_equal(cols,array([edge0,edge1,edge2,edge0,edge1,edge2,edge0,edge1,edge2]))


        rows,cols = massmatrix_rowcols(sc,2)

        assert_equal(rows,array([0]))
        assert_equal(cols,array([0]))


    
class test_whitney_innerproduct(TestCase):
    def setUp(self):
        self.cases = []

        #Self.Cases follow the pattern:
        #     (k,vertices,simplices,expected ordering of k-simplices,expected innerproduct matrix)        

        #0-form innerproduct for one unit line segment
        self.cases.append((0,array([[0],[1]]),array([[0,1]]),[(0,),(1,)], \
                           array([[1.0/3.0,1.0/6.0],[1.0/6.0,1.0/3.0]])))

        #0-form innerproduct for one unit line segment (reversed)
        self.cases.append((0,array([[1],[0]]),array([[1,0]]),[(1,),(0,)], \
                           array([[1.0/3.0,1.0/6.0],[1.0/6.0,1.0/3.0]])))

        #0-form innerproduct for one non-unit line segment
        self.cases.append((0,array([[-10],[3]]),array([[0,1]]),[(0,),(1,)], \
                           array([[13.0/3.0,13.0/6.0],[13.0/6.0,13.0/3.0]])))        

        #0-form innerproduct for two unit line segments
        self.cases.append((0,array([[0],[1],[2]]),array([[0,1],[1,2]]),[(0,),(1,),(2,)], \
                           array([[1.0/3.0,1.0/6.0,0.0],[1.0/6.0,2.0/3.0,1.0/6.0],[0.0,1.0/6.0,1.0/3.0]])))

        #0-form innerproduct for two non-unit line segments
        self.cases.append((0,array([[0],[3],[7.3]]),array([[0,1],[1,2]]),[(0,),(1,),(2,)], \
                           array([[3.0/3.0,3.0/6.0,0.0],[3.0/6.0,7.3/3.0,4.3/6.0],[0.0,4.3/6.0,4.3/3.0]])))

        #0-form innerproduct for reference triangle 
        self.cases.append((0,array([[0,0],[1,0],[0,1]]),array([[0,1,2]]),[(0,),(1,),(2,)], \
                            array([[1.0/12.0,1.0/24.0,1.0/24.0],[1.0/24.0,1.0/12.0,1.0/24.0],[1.0/24.0,1.0/24.0,1.0/12.0]]))) 

        #0-form innerproduct for reference tet
        self.cases.append((0,array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]),array([[0,1,2,3]]),[(0,),(1,),(2,),(3,)], \
                            array([[1.0/60.0,1.0/120.0,1.0/120.0,1.0/120.0], \
                                    [1.0/120.0,1.0/60.0,1.0/120.0,1.0/120.0], \
                                    [1.0/120.0,1.0/120.0,1.0/60.0,1.0/120.0], \
                                    [1.0/120.0,1.0/120.0,1.0/120.0,1.0/60.0]])))    

        #0-form innerproduct for reference tet (reversed)
        self.cases.append((0,array([[0,0,1],[0,1,0],[1,0,0],[0,0,0]]),array([[3,2,1,0]]),[(0,),(1,),(2,),(3,)], \
                            array([[1.0/60.0,1.0/120.0,1.0/120.0,1.0/120.0], \
                                    [1.0/120.0,1.0/60.0,1.0/120.0,1.0/120.0], \
                                    [1.0/120.0,1.0/120.0,1.0/60.0,1.0/120.0], \
                                    [1.0/120.0,1.0/120.0,1.0/120.0,1.0/60.0]])))    

        #1-form innerproduct for reference triangle
        self.cases.append((1,array([[0,0],[1,0],[0,1]]),array([[0,1,2]]),[(1,2),(0,1),(0,2)], \
                            array([[1.0/6.0,0.0,0.0],[0.0,1.0/3.0,1.0/6.0],[0.0,1.0/6.0,1.0/3.0]])))

        #1-form innerproduct for reference triangle (rotated vertices)
        self.cases.append((1,array([[1,0],[0,1],[0,0]]),array([[0,1,2]]),[(0,1),(0,2),(1,2)], \
                            array([[1.0/6.0,0.0,0.0],[0.0,1.0/3.0,1.0/6.0],[0.0,1.0/6.0,1.0/3.0]])))

        #1-form innerproduct for reference triangle (rotated indices)
        self.cases.append((1,array([[0,0],[1,0],[0,1]]),array([[1,2,0]]),[(1,2),(0,1),(0,2)], \
                            array([[1.0/6.0,0.0,0.0],[0.0,1.0/3.0,1.0/6.0],[0.0,1.0/6.0,1.0/3.0]])))


        #1-form innerproduct for scaled reference triangle
        self.cases.append((1,array([[0,0],[2,0],[0,2]]),array([[0,1,2]]),[(1,2),(0,1),(0,2)], \
                            array([[1.0/6.0,0.0,0.0],[0.0,1.0/3.0,1.0/6.0],[0.0,1.0/6.0,1.0/3.0]])))

        #1-form innerproduct for regular triangle
        self.cases.append((1,array([[0,0],[1,0],[0.5,sqrt(3.0)/2.0]]),array([[0,1,2]]),[(1,2),(0,1),(0,2)], \
                            array([[5*sqrt(3.0)/36.0, -sqrt(3.0)/36.0,  sqrt(3.0)/36.0], \
                                    [ -sqrt(3.0)/36.0,5*sqrt(3.0)/36.0,  sqrt(3.0)/36.0], \
                                    [  sqrt(3.0)/36.0,  sqrt(3.0)/36.0,5*sqrt(3.0)/36.0]])))

        #2-form innerproduct for reference triangle
        self.cases.append((2,array([[0,0],[1,0],[0,1]]),array([[0,1,2]]),[(0,1,2)], \
                            array([[12.0/6.0]])))

        #2-form innerproduct for reference triangle (rotated vertices)
        self.cases.append((2,array([[1,0],[0,1],[0,0]]),array([[0,1,2]]),[(0,1,2)], \
                            array([[12.0/6.0]])))

        #2-form innerproduct for reference triangle (rotated indices)
        self.cases.append((2,array([[0,0],[1,0],[0,1]]),array([[1,2,0]]),[(0,1,2)], \
                           array([[12.0/6.0]])))

        #1-form innerproduct for tetrahedron
        self.cases.append((1,array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]),array([[0,1,2,3]]),[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)], \
                            1.0/factorial(5)*array([[10,5,5,0,0,0], \
                                                     [5,10,5,0,0,0], \
                                                     [5,5,10,0,0,0], \
                                                     [0,0,0,4,1,-1], \
                                                     [0,0,0,1,4,1], \
                                                     [0,0,0,-1,1,4],])))

        #2-form innerproduct for reference tetrahedron
        self.cases.append((2,array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]),array([[0,1,2,3]]),[(1,2,3),(0,2,3),(0,1,3),(0,1,2)], \
                            4.0/factorial(5)*array([[6,-1, 1,-1], \
                                                    [-1,16, 4,-4], \
                                                    [ 1, 4,16,4], \
                                                    [-1,-4,4,16]])))

        #2-form innerproduct for reference tetrahedron (swapped vertices & indices)
        self.cases.append((2,array([[0,0,0],[0,1,0],[0,0,1],[1,0,0]]),array([[0,3,1,2]]),[(1,2,3),(0,1,2),(0,3,2),(0,3,1)], \
                            4.0/factorial(5)*array([[6, -1, -1,  1], \
                                                     [ -1, 16, -4,  4], \
                                                     [ -1, -4, 16,  4],
                                                     [  1,  4,  4, 16]])))

            
    def test_all(self):
        for k,v,s,o,expected in self.cases:
            sc = simplicial_complex((v,s))
            M = whitney_innerproduct(sc,k)
            M = M.todense()
            
            #Permute the matrix to coincide with the ordering of the known matrix
            permute = [sc[k].simplex_to_index[simplex(x)] for x in o]                                           

            M = M[permute,:][:,permute]
            
            assert_almost_equal(M,expected)
            
            #check whether matrix is S.P.D.
            self.assert_(alltrue(isreal(eigvals(M))))
            self.assert_(min(real(eigvals(M))) >= 0)
            assert_almost_equal(M,M.T)




class TestRegularCubeInnerproduct(TestCase):
    def setUp(self):
        pass


    def get_K(self,bitmap,k):
        """
        Receive a bitmap and return the k-form cube innerproduct
        """
        rcm = regular_cube_mesh(bitmap)
        rcc = regular_cube_complex(rcm)
        K = regular_cube_innerproduct(rcc,k).todense()
        return K
    

    def test_d1k0(self):
        #0-form innerproduct for unit line segment
        bitmap = ones((1,),dtype='bool')

        K_correct = array([[ 1.0/3.0,  1.0/6.0],
                           [ 1.0/6.0,  1.0/3.0]])

        assert_almost_equal(K_correct,self.get_K(bitmap,0))

    def test_d1k0(self):
        #1-form innerproduct for unit line segment
        bitmap = ones((1,),dtype='bool')

        K_correct = array([[ 1.00000000]])

        assert_almost_equal(K_correct,self.get_K(bitmap,1))


    def test_d2k0(self):
        #0-form innerproduct for unit square
        bitmap = ones((1,1),dtype='bool')


        K_correct = array([[  1.0/9.0,  1.0/18.0,  1.0/18.0,  1.0/36.0],
                           [ 1.0/18.0,   1.0/9.0,  1.0/36.0,  1.0/18.0],
                           [ 1.0/18.0,  1.0/36.0,   1.0/9.0,  1.0/18.0],
                           [ 1.0/36.0,  1.0/18.0,  1.0/18.0,   1.0/9.0]])

        assert_almost_equal(K_correct,self.get_K(bitmap,0))


    def test_d2k1(self):
        #1-form innerproduct for one unit square
        bitmap = ones((1,1),dtype='bool')

        K_correct = array([[ 1.0/3.0,        0,  1.0/6.0,        0],
                           [       0,  1.0/3.0,        0,  1.0/6.0],
                           [ 1.0/6.0,        0,  1.0/3.0,        0],
                           [       0,  1.0/6.0,        0,  1.0/3.0]])

        assert_almost_equal(K_correct,self.get_K(bitmap,1))


    def test_d2k1_2by2(self):
        #1-form innerproduct for 2x2 unit square grid
        bitmap = ones((2,2),dtype='bool')

        K_local = array([[ 1.0/3.0,        0,  1.0/6.0,        0],
                         [       0,  1.0/3.0,        0,  1.0/6.0],
                         [ 1.0/6.0,        0,  1.0/3.0,        0],
                         [       0,  1.0/6.0,        0,  1.0/3.0]])

        K_correct = zeros((12,12))

        K_correct[ix_([0,1,2,6],[0,1,2,6])]   = K_correct[ix_([0,1,2,6],[0,1,2,6])]   + K_local
        K_correct[ix_([5,6,7,10],[5,6,7,10])] = K_correct[ix_([5,6,7,10],[5,6,7,10])] + K_local
        K_correct[ix_([2,3,4,8],[2,3,4,8])]   = K_correct[ix_([2,3,4,8],[2,3,4,8])]   + K_local
        K_correct[ix_([7,8,9,11],[7,8,9,11])] = K_correct[ix_([7,8,9,11],[7,8,9,11])] + K_local

        assert_almost_equal(K_correct,self.get_K(bitmap,1))


    def test_d2k2(self):
        #1-form innerproduct for one unit square
        bitmap = ones((1,1),dtype='bool')

        K_correct = array([[ 1.0 ]])
        
        assert_almost_equal(K_correct,self.get_K(bitmap,2))


    def test_d3k0(self):
        #0-form innerproduct for unit cube
        bitmap = ones((1,1,1),dtype='bool')

        K_correct = array([[  1.0/27.0,  1.0/54.0,  1.0/54.0, 1.0/108.0,  1.0/54.0, 1.0/108.0, 1.0/108.0, 1.0/216.0],
                           [  1.0/54.0,  1.0/27.0, 1.0/108.0,  1.0/54.0, 1.0/108.0,  1.0/54.0, 1.0/216.0, 1.0/108.0],
                           [  1.0/54.0, 1.0/108.0,  1.0/27.0,  1.0/54.0, 1.0/108.0, 1.0/216.0,  1.0/54.0, 1.0/108.0],
                           [ 1.0/108.0,  1.0/54.0,  1.0/54.0,  1.0/27.0, 1.0/216.0, 1.0/108.0, 1.0/108.0,  1.0/54.0],
                           [  1.0/54.0, 1.0/108.0, 1.0/108.0, 1.0/216.0,  1.0/27.0,  1.0/54.0,  1.0/54.0, 1.0/108.0],
                           [ 1.0/108.0,  1.0/54.0, 1.0/216.0, 1.0/108.0,  1.0/54.0,  1.0/27.0, 1.0/108.0,  1.0/54.0],
                           [ 1.0/108.0, 1.0/216.0,  1.0/54.0, 1.0/108.0,  1.0/54.0, 1.0/108.0,  1.0/27.0,  1.0/54.0],
                           [ 1.0/216.0, 1.0/108.0, 1.0/108.0,  1.0/54.0, 1.0/108.0,  1.0/54.0,  1.0/54.0,  1.0/27.0]])

        assert_almost_equal(K_correct,self.get_K(bitmap,0))


    def test_d3k1(self):
        #1-form innerproduct for one unit square
        bitmap = ones((1,1,1),dtype='bool')

        K_correct = array([[  1.0/9.0,        0,        0, 1.0/18.0,        0, 1.0/18.0,        0, 1.0/36.0,        0,        0,        0,        0],
                           [        0,  1.0/9.0,        0,        0, 1.0/18.0,        0,        0,        0, 1.0/18.0,        0, 1.0/36.0,        0],
                           [        0,        0,  1.0/9.0,        0,        0,        0, 1.0/18.0,        0,        0, 1.0/18.0,        0, 1.0/36.0],
                           [ 1.0/18.0,        0,        0,  1.0/9.0,        0, 1.0/36.0,        0, 1.0/18.0,        0,        0,        0,        0],
                           [        0, 1.0/18.0,        0,        0,  1.0/9.0,        0,        0,        0, 1.0/36.0,        0, 1.0/18.0,        0],
                           [ 1.0/18.0,        0,        0, 1.0/36.0,        0,  1.0/9.0,        0, 1.0/18.0,        0,        0,        0,        0],
                           [        0,        0, 1.0/18.0,        0,        0,        0,  1.0/9.0,        0,        0, 1.0/36.0,        0, 1.0/18.0],
                           [ 1.0/36.0,        0,        0, 1.0/18.0,        0, 1.0/18.0,        0,  1.0/9.0,        0,        0,        0,        0],
                           [        0, 1.0/18.0,        0,        0, 1.0/36.0,        0,        0,        0,  1.0/9.0,        0, 1.0/18.0,        0],
                           [        0,        0, 1.0/18.0,        0,        0,        0, 1.0/36.0,        0,        0,  1.0/9.0,        0, 1.0/18.0],
                           [        0, 1.0/36.0,        0,        0, 1.0/18.0,        0,        0,        0, 1.0/18.0,        0,  1.0/9.0,        0],
                           [        0,        0, 1.0/36.0,        0,        0,        0, 1.0/18.0,        0,        0, 1.0/18.0,        0,  1.0/9.0]])
        
        assert_almost_equal(K_correct,self.get_K(bitmap,1))
        

    def test_d3k2(self):
        #1-form innerproduct for one unit square
        bitmap = ones((1,1,1),dtype='bool')

        K_correct = array([[ 1.0/3.0,       0,       0, 1.0/6.0,       0,       0],
                           [       0, 1.0/3.0,       0,       0, 1.0/6.0,       0],
                           [       0,       0, 1.0/3.0,       0,       0, 1.0/6.0],
                           [ 1.0/6.0,       0,       0, 1.0/3.0,       0,       0],
                           [       0, 1.0/6.0,       0,       0, 1.0/3.0,       0],
                           [       0,       0, 1.0/6.0,       0,       0, 1.0/3.0]])
        
        assert_almost_equal(K_correct,self.get_K(bitmap,2))

