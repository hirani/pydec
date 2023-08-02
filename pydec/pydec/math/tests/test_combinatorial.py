from pydec.testing import *

from scipy.special import factorial, comb

from pydec.math.combinatorial import combinations, permutations

def test_combinations():
    for N in xrange(6):
        L = range(N)            
        for K in xrange(N+1):                
            C = list(combinations(L,K))
            S = set([ frozenset(x) for x in C])
            ##Check order
            assert_equal(C, sorted(C))
            ##Check number of elements
            assert_equal(len(C), comb(N,K,exact=True))
            ##Make sure each element is unique
            assert_equal(len(S), comb(N,K,exact=True))
                   

def test_permutations():
    assert_equal(list(permutations([])),[])
    for N in xrange(1,6):
        L = range(N)            
        P = list(permutations(L))
        S = set([ frozenset(x) for x in P])
        ##Check number of elements
        assert_equal(len(P), factorial(N,exact=True))
        ##Check that there is only 1 unordered element
        assert_equal(len(S), 1)

