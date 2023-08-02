from pydec.testing import *

from pydec.math.parity import permutation_parity, relative_parity

cases = []
cases.append(([],0))  #for consistency
cases.append(([0],0))        
cases.append(([0,1],0))
cases.append(([1,0],1))
cases.append(([0,1,2],0))
cases.append(([0,2,1],1))
cases.append(([2,0,1],0))
cases.append(([2,1,0],1))
cases.append(([1,2,0],0))
cases.append(([1,0,2],1))
cases.append(([0,1,2,3],0))
cases.append(([0,1,3,2],1))
cases.append(([1,0,2,3],1))
cases.append(([1,2,0,3],0))
cases.append(([1,2,3,0],1))
cases.append(([2,1,3,0],0))
cases.append(([2,3,1,0],1))
cases.append(([3,2,1,0],0))
        
def test_permutation_parity():
    for perm,parity in cases:
        assert_equal(permutation_parity(perm), parity)

def test_relative_parity():
    for perm,parity in cases:
        #check lists of integers
        assert_equal(relative_parity(perm,perm), 0) 
        assert_equal(relative_parity(perm,range(len(perm))), parity)
        
        #check lists of strings        
        L1 = map(str,perm)
        L2 = map(str,range(len(perm)))
        assert_equal(relative_parity(L1,L1), 0)
        assert_equal(relative_parity(L1,L2), parity)

