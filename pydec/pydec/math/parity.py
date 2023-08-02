__all__ = ['relative_parity','permutation_parity']


def relative_parity(A,B):
    """Relative parity between two lists
    
    Parameters
    ----------
    A,B : lists of elements
        Lists A and B must contain permutations of the same elements.
    
        
    Returns
    -------
    parity : integer
        The parity is 0 if A differs from B by an even number of 
        transpositions and 1 otherwise.

    Examples
    --------
    >>> relative_parity( [0,1], [0,1] )
    0
    >>> relative_parity( [0,1], [1,0] )
    1
    >>> relative_parity( [0,1,2], [0,1,2] )
    0
    >>> relative_parity( [0,1,2], [0,2,1] )
    1
    >>> relative_parity( ['A','B','C'], ['A','B','C'] )
    0
    >>> relative_parity( ['A','B','C'], ['A','C','B'] )
    1

    """
    
    if len(A) != len(B): raise ValueError("B is not a permutation of A")
    
    # represent each element in B with its index in A and run permutation_parity()
    A_indices = dict(zip(A,range(len(A)))) 
    
    if len(A_indices) != len(A): raise ValueError("A contains duplicate values")
    
    try:   
        perm   = [A_indices[x] for x in B]
    except KeyError:
        raise ValueError("B is not a permutation of A")
    
    return permutation_parity(perm, check_input=False)
    
    

def permutation_parity(perm, check_input=True):
    """Parity of a permutation of the integers
    
    Parameters
    ----------
    perm : list of integers
        List containing a permutation of the integers 0...N
    
    Optional Parameters
    -------------------
    check_input : boolean
        If True, check whether the input is a valid permutation.
        
    Returns
    -------
    parity : integer
        The parity is 0 if perm differs from range(len(perm)) by an 
        even number of transpositions and 1 otherwise.

    Examples
    --------
    >>> permutation_parity( [0,1,2] )
    0
    >>> permutation_parity( [0,2,1] ) 
    1
    >>> permutation_parity( [1,0,2] )
    1
    >>> permutation_parity( [1,2,0] )
    0
    >>> permutation_parity( [2,0,1] )
    0
    >>> permutation_parity( [0,1,3,2] )
    1

    """

    n = len(perm)
    if check_input:
        rangen = range(n)
        if sorted(perm) != rangen:
            raise ValueError("Invalid input")

    # Decompose into disjoint cycles. We only need to
    # count the number of cycles to determine the parity
    num_cycles = 0
    seen = set()
    for i in range(n):
        if i in seen:
            continue
        num_cycles += 1
        j = i
        while True:
            assert j not in seen
            seen.add(j)
            j = perm[j]
            if j == i:
                break

    return (n - num_cycles) % 2

