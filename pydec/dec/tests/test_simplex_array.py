from pydec.testing import *

from scipy import random,arange,alltrue,array

from pydec.dec.simplex_array import simplex_array_boundary, \
     simplex_array_parity, simplex_array_searchsorted
from pydec.math.parity import relative_parity



class TestSearchsorted(TestCase):
    def setUp(self):
        random.seed(0)


    def test_simple1(self):
        s = array([[0],[1],[4],[6],[7],[10]])
        v = array([[6],[10],[0]])
        
        result   = simplex_array_searchsorted(s,v)
        expected = array([3,5,0])
        
        assert_equal(result,expected)


    def test_simple2(self):
        s = array([[0,1],[0,2],[1,2],[1,3],[1,4],[3,4]])
        v = array([[1,2],[0,2],[3,4]])
        
        result   = simplex_array_searchsorted(s,v)
        expected = array([2,1,5])

        assert_equal(result,expected)


    def test_random(self):
        for n_row in [1,2,3,10,100,200]:
            for n_col in [1,2,3,4,5]:
                s = arange(n_row*n_col).reshape((n_row,n_col))

                for n_searches in [1,2,3,n_row,2*n_row]:
                    expected = random.randint(0,n_row,n_searches)

                    v = s[expected,:]

                    result = simplex_array_searchsorted(s,v)

                    assert_equal(result,expected)

                


         
class test_simplex_array_parity(TestCase):
    def setUp(self):
        random.seed(0)

    def test_simple(self):
        cases = []
        cases.append(([[0]],[0]))
        cases.append(([[1]],[0]))
        cases.append(([[0],[1]],[0,0]))
        cases.append(([[0,1]],[0]))
        cases.append(([[1,0]],[1]))
        cases.append(([[0,1],[1,0]],[0,1]))
        cases.append(([[53,2]],[1]))
        cases.append(([[117,235]],[0]))
        cases.append(([[0,1,4]],[0]))
        cases.append(([[0,4,1]],[1]))
        cases.append(([[4,0,1]],[0]))
        cases.append(([[0,1,4],[0,4,1],[4,0,1]],[0,1,0]))

        for s,p in cases:
            assert_equal(simplex_array_parity(array(s)),array(p))

        
    def test_parity(self):
        """
        test parity with random data
        """
        for n_col in range(1,7):
            for n_row in range(1,50):
                A = arange(n_col*n_row)
                random.shuffle(A)
                A = A.reshape((n_row,n_col))

                s_parity = simplex_array_parity(A)

                for n,row in enumerate(A):
                    assert_equal(s_parity[n],relative_parity(row,sorted(row)))

                    
class TestBoundary(TestCase):

    def test_simple(self):
        cases = []
        cases.append(([[0,1]],[0],[[0],[1]],[[-1],[1]]))
        cases.append(([[0,1]],[1],[[0],[1]],[[1],[-1]]))
        cases.append(([[2,9]],[0],[[2],[9]],[[-1],[1]]))
        cases.append(([[2,9]],[1],[[2],[9]],[[1],[-1]]))
        cases.append(([[0,1,2]],[0],[[0,1],[0,2],[1,2]],[[1],[-1],[1]]))
        cases.append(([[0,1,2]],[1],[[0,1],[0,2],[1,2]],[[-1],[1],[-1]]))


        for simplex_array,parity,faces,boundary in cases:

            F,B = simplex_array_boundary(array(simplex_array),array(parity))

            assert_equal(F,array(faces))
            assert_equal(B.todense(),array(boundary))
                
