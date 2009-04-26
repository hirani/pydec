from pydec.testing import *

import numpy
from numpy import array, matrix
from scipy import rand
from scipy.linalg import norm

from pydec.dec.rips_complex import rips_complex, rips_simplices, \
        rips_chain_complex

def ensure_complex_exactness(cmplx):
    for d1,d2 in zip(cmplx.cochain_complex()[:-1], cmplx.cochain_complex()[1:]):
        assert_equal( (d2*d1).nnz, 0 )
    
    for b1,b2 in zip(cmplx.chain_complex()[:-1], cmplx.chain_complex()[1:]):
        assert_equal( (b1*b2).nnz, 0 )

    for d,b in zip(cmplx.cochain_complex()[:-1], cmplx.chain_complex()[1:]):
        assert_equal( (d.T - b).nnz, 0 )

    assert_equal( len(cmplx.chain_complex()), len(cmplx.cochain_complex()) )        


class TestRipsSimplices(TestCase):
    def setUp(self):
        numpy.random.seed(0)

    def test_simple1(self):
        """example with 1 triangle"""
        edges = array([[0,1],[0,2],[1,2]])

        expected = [ array([[0],[1],[2]]),
                     array([[0,1],[0,2],[1,2]]),
                     array([[0,1,2]]) ]

        result = rips_simplices(3, edges, 2)

        for r,e in zip(result,expected):
            assert_equal(r,e)

    def test_simple2(self):
        """example with 1 tet and 1 triangle"""
        edges = array([[0,1],[0,2],[0,3],[0,4],[0,5],[1,2],[1,5],[2,5],[3,4]])

        expected = [ array([[0],[1],[2],[3],[4],[5]]),
                     array([[0,1],[0,2],[0,3],
                            [0,4],[0,5],[1,2],
                            [1,5],[2,5],[3,4]]),
                     array([[0,1,2],[0,1,5],[0,2,5],[0,3,4],[1,2,5]]),
                     array([[0,1,2,5]]) ]

        result = rips_simplices(6, edges, 3)

        for r,e in zip(result,expected):
            assert_equal(r,e)

    def test_random_2d(self):
        N = 200
        R = 0.2
    
        pts = rand(N,2)
        rc  = rips_complex( pts, R )
    
        edges = set( [tuple(e) for e in rc.simplices[1]] )
        
        for i in range(N):
            for j in range(i+1,N):
                if norm(pts[i] - pts[j]) < R:
                    assert( (i,j) in edges )
                else:
                    assert( (i,j) not in edges )
    
        for t in rc.simplices[2]:
            for i,j in [(0,1),(0,2),(1,2)]:
                assert( (t[i],t[j]) in edges )
    
        ensure_complex_exactness(rc)
        
    def test_random_3d(self):
        N = 200
        R = 0.25
    
        pts = rand(N,3)
        rc  = rips_complex( pts, R )
    
        edges = set( [tuple(e) for e in rc.simplices[1]] )
        
        for i in range(N):
            for j in range(i+1,N):
                if norm(pts[i] - pts[j]) < R:
                    assert( (i,j) in edges )
                else:
                    assert( (i,j) not in edges )
       
        for t in rc.simplices[2]:
            for i,j in [(0,1),(0,2),(1,2)]:
                assert( (t[i],t[j]) in edges )
    
        for t in rc.simplices[3]:
            for i,j in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]:
                assert( (t[i],t[j]) in edges )
        
        ensure_complex_exactness(rc)


class TestRipsComplex(TestCase):
    def test_simple1(self):
        """example with 1 edge and 1 point"""

        simplices = [ array([[0],[1],[2]]),
                      array([[0,2]]) ]

        expected = [ matrix([[0,0,0]]),
                     matrix([[-1],[ 0],[ 1]]) ]

        result = rips_chain_complex(simplices)

        for r,e in zip(result,expected):
            assert_equal(r.todense(),e)

    def test_simple2(self):
        """example with 1 triangle"""

        simplices = [ array([[0],[1],[2]]),
                      array([[0,1],[0,2],[1,2]]),
                      array([[0,1,2]]) ]

        expected = [ matrix([[0,0,0]]),
                     matrix([[-1,-1, 0],
                             [ 1, 0,-1],
                             [ 0, 1, 1]]),
                     matrix([[ 1],[-1],[ 1]]) ]

        result = rips_chain_complex(simplices)

        for r,e in zip(result,expected):
            assert_equal(r.todense(),e)

    def test_simple3(self):
        """example with 2 triangles and 1 edge"""

        simplices = [ array([[0],[1],[2],[3],[4]]),
                      array([[0, 1],
                             [0, 2],
                             [0, 3],
                             [1, 2],
                             [2, 3],
                             [2, 4]]),
                      array([[0, 1, 2],
                             [0, 2, 3]]) ]

        expected = [ array([[ 0, 0, 0, 0, 0]]),
                     array([[-1,-1,-1, 0, 0, 0],
                            [ 1, 0, 0,-1, 0, 0],
                            [ 0, 1, 0, 1,-1,-1],
                            [ 0, 0, 1, 0, 1, 0],
                            [ 0, 0, 0, 0, 0, 1]]),
                     array([[ 1, 0],
                            [-1, 1],
                            [ 0,-1],
                            [ 1, 0],
                            [ 0, 1],
                            [ 0, 0]]) ]
        
        result = rips_chain_complex( simplices )
        
        for r,e in zip(result,expected):
            assert_equal(r.todense(),e)
                     


    def test_simple4(self):
        """example with 1 triangle and 1 edge"""

        simplices = [ array([[0],[1],[2],[3]]),
                      array([[0,1],[0,2],[0,3],[1,2]]),
                      array([[0,1,2]]) ]

        expected = [ matrix([[ 0, 0, 0, 0]]),
                     matrix([[-1,-1,-1, 0],
                             [ 1, 0, 0,-1],
                             [ 0, 1, 0, 1],
                             [ 0, 0, 1, 0]]),
                     matrix([[ 1],[-1],[ 0],[ 1]]) ]

        result = rips_chain_complex( simplices )

        for r,e in zip(result,expected):
            assert_equal(r.todense(),e)

