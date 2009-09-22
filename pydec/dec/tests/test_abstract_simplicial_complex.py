from pydec.testing import *

from pydec.dec import abstract_simplicial_complex
         
class TestAbstractSimplicialComplex(TestCase):
    def test_1d(self):
        # two edges
        s = [  [[0,1],[2,3]]  ]
        asc = abstract_simplicial_complex(s)
        assert_equal(asc.simplices[0], [[0],[1],[2],[3]])
        assert_equal(asc.simplices[1], [[0,1],[2,3]])
        
        B0 = asc.chain_complex()[0].todense()
        B1 = asc.chain_complex()[1].todense()

        assert_equal(B0, [[0,0,0,0]])
        assert_equal(B1, [[-1, 0],
                          [ 1, 0],
                          [ 0,-1],
                          [ 0, 1]])
        
        # two edges, two redundant vertices
        s = [  [[0],[2]],  [[0,1],[2,3]]  ]
        asc = abstract_simplicial_complex(s)
        assert_equal(asc.simplices[0], [[0],[1],[2],[3]])
        assert_equal(asc.simplices[1], [[0,1],[2,3]])
        
        B0 = asc.chain_complex()[0].todense()
        B1 = asc.chain_complex()[1].todense()

        assert_equal(B0, [[0,0,0,0]])
        assert_equal(B1, [[-1, 0],
                          [ 1, 0],
                          [ 0,-1],
                          [ 0, 1]])

        # two edges, two maximal vertices
        s = [  [[2],[5]],  [[0,1],[3,4]]  ]
        asc = abstract_simplicial_complex(s)
        assert_equal(asc.simplices[0], [[0],[1],[2],[3],[4],[5]])
        assert_equal(asc.simplices[1], [[0,1],[3,4]])
        
        B0 = asc.chain_complex()[0].todense()
        B1 = asc.chain_complex()[1].todense()

        assert_equal(B0, [[0,0,0,0,0,0]])
        assert_equal(B1, [[-1, 0],
                          [ 1, 0],
                          [ 0, 0],
                          [ 0,-1],
                          [ 0, 1],
                          [ 0, 0]])
    
    def test_2d(self):
        # one triangle
        s = [  [[0,1,2]]  ]
        asc = abstract_simplicial_complex(s)
        assert_equal(asc.simplices[0], [[0],[1],[2]])
        assert_equal(asc.simplices[1], [[0,1],[0,2],[1,2]])
        assert_equal(asc.simplices[2], [[0,1,2]])
        
        B0 = asc.chain_complex()[0].todense()
        B1 = asc.chain_complex()[1].todense()
        B2 = asc.chain_complex()[2].todense()

        assert_equal(B0, [[0,0,0]])
        assert_equal(B1, [[-1,-1, 0],
                          [ 1, 0,-1],
                          [ 0, 1, 1]])
        assert_equal(B2, [[ 1],
                          [-1],
                          [ 1]])
        
        # one triangle, one redundant edge
        s = [  [[0,2]],  [[0,1,2]]  ]
        asc = abstract_simplicial_complex(s)
        assert_equal(asc.simplices[0], [[0],[1],[2]])
        assert_equal(asc.simplices[1], [[0,1],[0,2],[1,2]])
        assert_equal(asc.simplices[2], [[0,1,2]])
        
        B0 = asc.chain_complex()[0].todense()
        B1 = asc.chain_complex()[1].todense()
        B2 = asc.chain_complex()[2].todense()

        assert_equal(B0, [[0,0,0]])
        assert_equal(B1, [[-1,-1, 0],
                          [ 1, 0,-1],
                          [ 0, 1, 1]])
        assert_equal(B2, [[ 1],
                          [-1],
                          [ 1]])
        
        # one triangle, two maximal edges
        s = [  [[0,3],[4,5]],  [[0,1,2]]  ]
        asc = abstract_simplicial_complex(s)
        assert_equal(asc.simplices[0], [[0],[1],[2],[3],[4],[5]])
        assert_equal(asc.simplices[1], [[0,1],[0,2],[0,3],[1,2],[4,5]])
        assert_equal(asc.simplices[2], [[0,1,2]])
        
        B0 = asc.chain_complex()[0].todense()
        B1 = asc.chain_complex()[1].todense()
        B2 = asc.chain_complex()[2].todense()

        assert_equal(B0, [[0,0,0,0,0,0]])
        assert_equal(B1, [[-1,-1,-1, 0, 0],
                          [ 1, 0, 0,-1, 0],
                          [ 0, 1, 0, 1, 0],
                          [ 0, 0, 1, 0, 0],
                          [ 0, 0, 0, 0,-1],
                          [ 0, 0, 0, 0, 1]])
        assert_equal(B2, [[ 1],
                          [-1],
                          [ 0],
                          [ 1],
                          [ 0]])


