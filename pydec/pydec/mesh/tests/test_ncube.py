from pydec.testing import *

from numpy import arange, prod, array

from pydec.mesh.ncube import nCube

  
class test_relative_partiy(TestCase):
    def setUp(self):
        test_cases = []
        test_cases.append(([0],0))
        test_cases.append(([0,1],0))
        test_cases.append(([1,0],1))
        test_cases.append(([[0,1],[2,3]],0))
        test_cases.append(([[1,3],[0,2]],0))
        test_cases.append(([[3,2],[1,0]],0))
        test_cases.append(([[2,0],[3,1]],0))
        test_cases.append(([[2,3],[0,1]],1))
        test_cases.append(([[0,2],[1,3]],1))
        test_cases.append(([[1,0],[3,2]],1))
        test_cases.append(([[3,1],[2,0]],1))
        test_cases.append(([[3,1],[2,0]],1))
        test_cases.append(([[[0,1],[2,3]],[[4,5],[6,7]]],0))
        test_cases.append(([[[1,5],[3,7]],[[0,4],[2,6]]],0))
        test_cases.append(([[[5,4],[7,6]],[[1,0],[3,2]]],0))
        test_cases.append(([[[4,0],[6,2]],[[5,1],[7,3]]],0))
        test_cases.append(([[[2,3],[6,7]],[[0,1],[4,5]]],0))
        test_cases.append(([[[4,5],[0,1]],[[6,7],[2,3]]],0))
        test_cases.append(([[[7,3],[5,1]],[[6,2],[4,0]]],0))
        test_cases.append(([[[1,0],[5,4]],[[3,2],[7,6]]],0))
        test_cases.append(([[[4,5],[6,7]],[[0,1],[2,3]]],1))
        test_cases.append(([[[0,4],[2,6]],[[1,5],[3,7]]],1))
        test_cases.append(([[[1,0],[3,2]],[[5,4],[7,6]]],1))
        test_cases.append(([[[5,1],[7,3]],[[4,0],[6,2]]],1))
        self.test_cases = test_cases
        

    def check_basic(self):
        """
        Test permutations relative to the canonical order
        """
        for B,result in self.test_cases:
            B = array(B)
            A = arange(prod(B.shape)).reshape(B.shape)
            
            A_cube = nCube(A)
            B_cube = nCube(B)

            assert_equal(A_cube.relative_orientation(A_cube),False)
            assert_equal(B_cube.relative_orientation(B_cube),False)
            assert_equal(A_cube.relative_orientation(B_cube),result)
            assert_equal(B_cube.relative_orientation(A_cube),result)

    def check_pairs(self):
        """
        Test pairs of cubes against one another
        """
        
        for B0,result0 in self.test_cases:
            B0 = array(B0)
            for B1,result1 in self.test_cases:
                B1 = array(B1)
                if B0.shape == B1.shape:
                    B0_cube = nCube(B0)
                    B1_cube = nCube(B1)
                    
                    assert_equal(B0_cube.relative_orientation(B1_cube),result0 ^ result1)
                    assert_equal(B1_cube.relative_orientation(B0_cube),result0 ^ result1)


    def check_permutations(self):
        """
        Test permuted versions of the canonical ordering
        """
        for N in range(2,6):
            A = arange(2**N).reshape((2,)*N)
            B = A.copy()
            for i in range(1,N+1):
                Bcopy = B.copy()
                slice0 = [slice(None,None,None)]*(i-1) + [0] + [slice(None,None,None)]*(N-i)
                slice1 = [slice(None,None,None)]*(i-1) + [1] + [slice(None,None,None)]*(N-i)
                B[slice0] = Bcopy[slice1]
                B[slice1] = Bcopy[slice0]

                A_cube = nCube(A)
                B_cube = nCube(B)
                assert_equal(A_cube.relative_orientation(B_cube),i % 2)
                assert_equal(B_cube.relative_orientation(A_cube),i % 2)

        
