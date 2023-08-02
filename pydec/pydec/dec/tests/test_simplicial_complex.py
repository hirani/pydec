from pydec.testing import *

from scipy import random, rand, concatenate, zeros, sparse, matrix, sqrt

from pydec.dec import simplicial_complex
from pydec.mesh.simplex import simplex

         
class TestSimplicialComplex(TestCase):
    def test_d_in_1d(self):
        #sorted test
        v,e = matrix([[0],[1]]),matrix([[0,1]])
        sc  = simplicial_complex((v,e))
        assert_equal(sc[0].d.shape,(1,2))
        assert_equal(sc[0].d[0,0],-1)
        assert_equal(sc[0].d[0,1], 1)

        #reversed test
        v,e = matrix([[0],[1]]),matrix([[1,0]])
        sc  = simplicial_complex((v,e))
        assert_equal(sc[0].d.shape,(1,2))
        assert_equal(sc[0].d[0,0], 1)
        assert_equal(sc[0].d[0,1],-1)


    def test_d_one_triangle(self):
        #sorted test
        v,e = matrix([[0,0],[1,0],[0,1]]),matrix([[0,1,2]])
        sc  = simplicial_complex((v,e))
        assert_equal(sc[0].d.shape,(3,3))

        #d0
        row = sc[1].simplex_to_index[simplex((0,1))]
        assert_equal(sc[0].d[row,0],-1)
        assert_equal(sc[0].d[row,1], 1)

        row = sc[1].simplex_to_index[simplex((1,2))]
        assert_equal(sc[0].d[row,1],-1)
        assert_equal(sc[0].d[row,2], 1)

        row = sc[1].simplex_to_index[simplex((0,2))]
        assert_equal(sc[0].d[row,0],-1)
        assert_equal(sc[0].d[row,2], 1)

        #d1
        col = sc[1].simplex_to_index[simplex((0,1))]
        assert_equal(sc[1].d[0,col], 1)

        col = sc[1].simplex_to_index[simplex((0,2))]
        assert_equal(sc[1].d[0,col],-1)

        col = sc[1].simplex_to_index[simplex((1,2))]
        assert_equal(sc[1].d[0,col], 1)

        #reversed test
        v,e = matrix([[0,0],[1,0],[0,1]]),matrix([[0,2,1]])
        sc  = simplicial_complex((v,e))
        assert_equal(sc[0].d.shape,(3,3))

        #d0
        row = sc[1].simplex_to_index[simplex((0,1))]
        assert_equal(sc[0].d[row,0],-1)
        assert_equal(sc[0].d[row,1], 1)

        row = sc[1].simplex_to_index[simplex((1,2))]
        assert_equal(sc[0].d[row,1],-1)
        assert_equal(sc[0].d[row,2], 1)

        row = sc[1].simplex_to_index[simplex((0,2))]
        assert_equal(sc[0].d[row,0],-1)
        assert_equal(sc[0].d[row,2], 1)

        #d1
        col = sc[1].simplex_to_index[simplex((0,1))]
        assert_equal(sc[1].d[0,col],-1)

        col = sc[1].simplex_to_index[simplex((0,2))]
        assert_equal(sc[1].d[0,col], 1)

        col = sc[1].simplex_to_index[simplex((1,2))]
        assert_equal(sc[1].d[0,col],-1)



    def test_d_two_triangles(self):
        v,e = matrix([[0,0],[1,0],[0,1],[1,1]]),matrix([[1,3,2],[2,0,1]])
        sc  = simplicial_complex((v,e))
        
        assert_equal(sc[0].d.shape,(5,4))
        assert_equal(sc[1].d.shape,(2,5))       

        #d0
        edge01 = sc[1].simplex_to_index[simplex((0,1))]
        edge02 = sc[1].simplex_to_index[simplex((0,2))]
        edge12 = sc[1].simplex_to_index[simplex((1,2))]
        edge13 = sc[1].simplex_to_index[simplex((1,3))]
        edge23 = sc[1].simplex_to_index[simplex((2,3))]        
        
        assert_equal(sc[0].d[edge01,0],-1)
        assert_equal(sc[0].d[edge01,1], 1)

        assert_equal(sc[0].d[edge02,0],-1)
        assert_equal(sc[0].d[edge02,2], 1)

        assert_equal(sc[0].d[edge12,1],-1)
        assert_equal(sc[0].d[edge12,2], 1)

        assert_equal(sc[0].d[edge13,1],-1)
        assert_equal(sc[0].d[edge13,3], 1)

        assert_equal(sc[0].d[edge23,2],-1)
        assert_equal(sc[0].d[edge23,3], 1)

        assert_equal(sc[1].d[1,edge01], 1)
        assert_equal(sc[1].d[1,edge12], 1)
        assert_equal(sc[1].d[1,edge02],-1)

        assert_equal(sc[1].d[0,edge23],-1)
        assert_equal(sc[1].d[0,edge12],-1)
        assert_equal(sc[1].d[0,edge13], 1)


    def test_d_one_tet(self):
        v,e = matrix([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]),matrix([[3,2,1,0]])
        sc  = simplicial_complex((v,e))
        
        assert_equal(sc[0].d.shape,(6,4))
        assert_equal(sc[1].d.shape,(4,6))
        assert_equal(sc[2].d.shape,(1,4))

        #d0
        edge01 = sc[1].simplex_to_index[simplex((0,1))]
        edge02 = sc[1].simplex_to_index[simplex((0,2))]
        edge03 = sc[1].simplex_to_index[simplex((0,3))]
        edge12 = sc[1].simplex_to_index[simplex((1,2))]
        edge13 = sc[1].simplex_to_index[simplex((1,3))]
        edge23 = sc[1].simplex_to_index[simplex((2,3))]        
        
        assert_equal(sc[0].d[edge01,0],-1)
        assert_equal(sc[0].d[edge01,1], 1)

        assert_equal(sc[0].d[edge02,0],-1)
        assert_equal(sc[0].d[edge02,2], 1)

        assert_equal(sc[0].d[edge03,0],-1)
        assert_equal(sc[0].d[edge03,3], 1)

        assert_equal(sc[0].d[edge12,1],-1)
        assert_equal(sc[0].d[edge12,2], 1)

        assert_equal(sc[0].d[edge13,1],-1)
        assert_equal(sc[0].d[edge13,3], 1)

        assert_equal(sc[0].d[edge23,2],-1)
        assert_equal(sc[0].d[edge23,3], 1)

        face123 = sc[2].simplex_to_index[simplex((1,2,3))]
        face023 = sc[2].simplex_to_index[simplex((0,2,3))]
        face013 = sc[2].simplex_to_index[simplex((0,1,3))]
        face012 = sc[2].simplex_to_index[simplex((0,1,2))]

        #d1
        assert_equal(sc[1].d[face012,edge12], 1)
        assert_equal(sc[1].d[face012,edge02],-1)
        assert_equal(sc[1].d[face012,edge01], 1)

        assert_equal(sc[1].d[face023,edge23], 1)
        assert_equal(sc[1].d[face023,edge03],-1)
        assert_equal(sc[1].d[face023,edge02], 1)

        assert_equal(sc[1].d[face013,edge13], 1)
        assert_equal(sc[1].d[face013,edge03],-1)
        assert_equal(sc[1].d[face013,edge01], 1)

        assert_equal(sc[1].d[face012,edge12], 1)
        assert_equal(sc[1].d[face012,edge02],-1)
        assert_equal(sc[1].d[face012,edge12], 1)

        #d2
        assert_equal(sc[2].d[0,face123], 1)
        assert_equal(sc[2].d[0,face023],-1)
        assert_equal(sc[2].d[0,face013], 1)
        assert_equal(sc[2].d[0,face012],-1)


    def test_d_two_tets(self):
        """
        Check d0,d1,d2 for a mesh with two tets
        """
        v,e = matrix([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1]]),matrix([[0,1,2,3],[1,2,3,4]])
        sc  = simplicial_complex((v,e))
        
        assert_equal(sc[0].d.shape,(9,5))
        assert_equal(sc[1].d.shape,(7,9))
        assert_equal(sc[2].d.shape,(2,7))

        #d0
        edge01 = sc[1].simplex_to_index[simplex((0,1))]
        edge02 = sc[1].simplex_to_index[simplex((0,2))]
        edge03 = sc[1].simplex_to_index[simplex((0,3))]
        edge12 = sc[1].simplex_to_index[simplex((1,2))]
        edge13 = sc[1].simplex_to_index[simplex((1,3))]
        edge23 = sc[1].simplex_to_index[simplex((2,3))]        
        edge14 = sc[1].simplex_to_index[simplex((1,4))]
        edge24 = sc[1].simplex_to_index[simplex((2,4))]
        edge34 = sc[1].simplex_to_index[simplex((3,4))]
        
        assert_equal(sc[0].d[edge01,0],-1)
        assert_equal(sc[0].d[edge01,1], 1)

        assert_equal(sc[0].d[edge02,0],-1)
        assert_equal(sc[0].d[edge02,2], 1)

        assert_equal(sc[0].d[edge03,0],-1)
        assert_equal(sc[0].d[edge03,3], 1)

        assert_equal(sc[0].d[edge12,1],-1)
        assert_equal(sc[0].d[edge12,2], 1)

        assert_equal(sc[0].d[edge13,1],-1)
        assert_equal(sc[0].d[edge13,3], 1)

        assert_equal(sc[0].d[edge23,2],-1)
        assert_equal(sc[0].d[edge23,3], 1)

        assert_equal(sc[0].d[edge14,1],-1)
        assert_equal(sc[0].d[edge14,4], 1)

        assert_equal(sc[0].d[edge24,2],-1)
        assert_equal(sc[0].d[edge24,4], 1)

        assert_equal(sc[0].d[edge34,3],-1)
        assert_equal(sc[0].d[edge34,4], 1)

        face123 = sc[2].simplex_to_index[simplex((1,2,3))]
        face023 = sc[2].simplex_to_index[simplex((0,2,3))]
        face013 = sc[2].simplex_to_index[simplex((0,1,3))]
        face012 = sc[2].simplex_to_index[simplex((0,1,2))]
        face124 = sc[2].simplex_to_index[simplex((1,2,4))]
        face134 = sc[2].simplex_to_index[simplex((1,3,4))]
        face234 = sc[2].simplex_to_index[simplex((2,3,4))]

        #d1
        assert_equal(sc[1].d[face012,edge12], 1)
        assert_equal(sc[1].d[face012,edge02],-1)
        assert_equal(sc[1].d[face012,edge01], 1)

        assert_equal(sc[1].d[face023,edge23], 1)
        assert_equal(sc[1].d[face023,edge03],-1)
        assert_equal(sc[1].d[face023,edge02], 1)

        assert_equal(sc[1].d[face013,edge13], 1)
        assert_equal(sc[1].d[face013,edge03],-1)
        assert_equal(sc[1].d[face013,edge01], 1)

        assert_equal(sc[1].d[face012,edge12], 1)
        assert_equal(sc[1].d[face012,edge02],-1)
        assert_equal(sc[1].d[face012,edge12], 1)

        assert_equal(sc[1].d[face124,edge24], 1)
        assert_equal(sc[1].d[face124,edge14],-1)
        assert_equal(sc[1].d[face124,edge12], 1)

        assert_equal(sc[1].d[face134,edge34], 1)
        assert_equal(sc[1].d[face134,edge14],-1)
        assert_equal(sc[1].d[face134,edge13], 1)

        assert_equal(sc[1].d[face234,edge34], 1)
        assert_equal(sc[1].d[face234,edge24],-1)
        assert_equal(sc[1].d[face234,edge23], 1)

        #d2
        assert_equal(sc[2].d[0,face123], 1)
        assert_equal(sc[2].d[0,face023],-1)
        assert_equal(sc[2].d[0,face013], 1)
        assert_equal(sc[2].d[0,face012],-1)

        assert_equal(sc[2].d[1,face234], 1)
        assert_equal(sc[2].d[1,face134],-1)
        assert_equal(sc[2].d[1,face124], 1)
        assert_equal(sc[2].d[1,face123],-1)

                               
    def test_hodge_star(self):
        """Test the hodge * operator"""
        regular_cases = []
        #Unit line segment
        regular_cases.append((matrix([[0],[1]]),matrix([[0,1]]),[1,1],[0.5,1]))
        #One regular triangle, edge length 1
        regular_cases.append((matrix([[0,0],[1,0],[0.5,sqrt(3.0)/2.0]]),matrix([[0,1,2]]), [1, 1, sqrt(3.0)/4.0],\
                            [sqrt(3.0)/12.0,sqrt(3.0)/6.0,1]))
        # One regular tet, edge length sqrt(2)
        regular_cases.append((matrix([[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 0]]),matrix([[0,1,2,3]]),\
                            [1, sqrt(2.0),sqrt(3.0/4.0), 1.0/3.0],[1.0/12.0,sqrt(2.0)/12.0,sqrt(3.0/36.0),1]))
        
        for v,e,primal,dual in regular_cases:
            sc = simplicial_complex((v,e))
            for dim in range(sc.complex_dimension() + 1):
                (n,m) = sc[dim].star.shape
                for i in range(n):                    
                    assert_almost_equal( sc[dim].star[i,i], dual[dim]/primal[dim] )
        
        
        multi_cases = []
        #two line segments
        multi_cases.append((matrix([[0],[1],[2]]),matrix([[0,1],[1,2]]), \
                        [([0],0.5),([1],1.0),([2],0.5),([0,1],1),([1,2],1)]))
        #three line segments
        multi_cases.append((matrix([[0],[1],[2],[3]]),matrix([[0,1],[1,2],[2,3]]), \
                        [([0],0.5),([1],1.0),([2],1),([3],0.5),([0,1],1),([1,2],1),([2,3],1)]))
        #two triangles
        multi_cases.append((matrix([[0,0],[1.0,0],[0.5,sqrt(3.0)/2.0],[1.5,sqrt(3.0)/2.0]]), \
                        matrix([[0,1,2],[1,3,2]]), \
                        [([0],sqrt(3.0)/12.0),([1],2*sqrt(3.0)/12.0),([2],2*sqrt(3.0)/12.0), ([3],sqrt(3.0)/12.0), \
                        ([0,1],sqrt(3.0)/6.0),([0,2],sqrt(3.0)/6.0),([1,2],2*sqrt(3.0)/6.0),([1,3],sqrt(3.0)/6.0), \
                        ([2,3],sqrt(3.0)/6.0),([0,1,2],4.0/sqrt(3.0)),([1,2,3],4.0/sqrt(3.0))]))
        #three triangles
        multi_cases.append((matrix([[0,0],[1.0,0],[0.5,sqrt(3.0)/2.0],[1.5,sqrt(3.0)/2.0],[2,0]]), \
                        matrix([[0,1,2],[1,3,2],[1,4,3]]), \
                        [([0],sqrt(3.0)/12.0),([1],3*sqrt(3.0)/12.0),([2],2*sqrt(3.0)/12.0), ([3],2*sqrt(3.0)/12.0), \
                        ([4],sqrt(3.0)/12.0),([0,1],sqrt(3.0)/6.0),([0,2],sqrt(3.0)/6.0),([1,2],2*sqrt(3.0)/6.0),
                        ([1,3],2*sqrt(3.0)/6.0), ([1,4],sqrt(3.0)/6.0),([3,4],sqrt(3.0)/6.0),\
                        ([2,3],sqrt(3.0)/6.0),([0,1,2],4.0/sqrt(3.0)),([1,2,3],4.0/sqrt(3.0)),([1,4,3],4.0/sqrt(3.0))]))
        
        
        for v,e,t_list in multi_cases:        
            sc = simplicial_complex((v,e))
            for s,ratio in t_list:
                data = sc[len(s)-1]
                index = data.simplex_to_index[simplex(s)]
                assert_almost_equal(ratio,data.star[index,index])
              
            
        #use the same multi_cases, but add dimensions and translate the vertices    
        #TODO rotate also
        random.seed(0)
        for v,e,t_list in multi_cases:
        
            (rows,cols) = v.shape
            v = concatenate([v,zeros((rows,4))],1)            
            v += rand(1,cols +4) #translate
            
            sc = simplicial_complex((v,e))
            for s,ratio in t_list:
                data = sc[len(s)-1]
                index = data.simplex_to_index[simplex(s)]
                assert_almost_equal(ratio,data.star[index,index])
            
        #Test sign of dual star,  check that:   ** = -1 ^(k (n-k))
        for v,e,t_list in multi_cases:
            sc = simplicial_complex((v,e))
            for dim,data in enumerate(sc):
                n = sc.complex_dimension()
                k = dim
                assert_almost_equal((data.star * data.star_inv).todense(), \
                                    ((-1)**(k*(n-k))*sparse.identity(data.num_simplices)).todense())
    
