__all__ = ['combinations','permutations']


def combinations(L, n):
    """
    Returns a generator object that iterates over 
    all combinations of n elements from a list L
    
    The elements remain in the same order as in L    
    """
    if n==0 or n > len(L): 
        yield []
    else:
        for i in xrange(len(L)-n+1):
            for t in combinations(L[i+1:],n-1):
                yield [L[i]]+t
        
def permutations(L):
    """
    Returns a generator object that iterates over 
    all permutations of elements in a list L
    """
    if len(L) == 1:
        yield [L[0]]
    elif len(L) >= 2:
        (h, t) = (L[0:1], L[1:])
        for p in permutations(t):
            for i in range(len(p)+1):
                yield p[:i] + h + p[i:]
 
 
          
