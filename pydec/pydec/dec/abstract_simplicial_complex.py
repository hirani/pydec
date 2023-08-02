from numpy import zeros, ones, arange, array, asarray, hstack, vstack, empty, lexsort, atleast_2d, alltrue
from scipy import sparse

from .simplex_array import simplex_array_parity, simplex_array_boundary, simplex_array_searchsorted

__all__ = ['abstract_simplicial_complex']

class abstract_simplicial_complex:
    def __init__(self, simplices):
        """Construct an abstract simplicial complex

        Parameters
        ----------
        simplices : list of arrays
            Maximal simplices of each dimension
            TODO

        Examples
        --------
        >>> from pydec.dec import abstract_simplicial_complex
        >>> from numpy import array
        >>> simplices = [array([[4]]), array([[0,3]])]
        >>> asc = abstract_simplicial_complex(simplices)

        TODO

        >>> print rc.simplices[0]
        >>> print rc.simplices[1]

        Notes
        -----

        TODO explain input handling

        """

        # convert array-like objects to arrays
        simplices = [atleast_2d(s) for s in simplices]

        # find top simplex dimension
        D = max([s.shape[1] for s in simplices]) - 1
         
        # convert simplices list to canonical simplex_array list representation
        old_simplices = simplices
        simplices = [None] * (D + 1)
        for s in old_simplices:
            simplices[s.shape[1] - 1] = s
       

        # top most simplex array
        s = simplices[-1].copy()

        parity = simplex_array_parity(s)
        s.sort()

        chain_complex = [None] * (D + 1)

        for d in range(D - 1, -1, -1):
            s,B = simplex_array_boundary(s,parity)

            if simplices[d] is not None:
                old_s = s

                # sort columns to ensure user-defined faces are in canonical format
                simplices[d].sort()

                # merge user-defined faces with boundary faces
                s = vstack((s,simplices[d]))
                
                # sort rows to bring equivalent elements together
                s = s[lexsort(s.T[::-1])]

                # find unique simplices
                mask = ~hstack((array([False]),alltrue(s[1:] == s[:-1],axis=1)))
                s = s[mask]

                # indices of the boundary faces in the full face array
                remap = simplex_array_searchsorted(s, old_s)
                
                # remap indices of boundary operator
                B = B.tocoo(copy=False)
                B = sparse.coo_matrix((B.data, (remap[B.row], B.col)), (s.shape[0], B.shape[1]))
                B = B.tocsr()
            
            # faces are already in canonical format, so parity is even
            parity = zeros(s.shape[0],dtype=s.dtype)
            
            simplices[d]       = s
            chain_complex[d+1] = B

        # compute 0-simplices and boundary operator
        simplices[0]     = arange(simplices[0].max() + 1).reshape(-1,1)
        chain_complex[0] = sparse.csr_matrix( (1,len(s)), dtype='uint8')

        # store the cochain complex
        Bn = chain_complex[-1]
        cochain_complex  = [ B.T for B in chain_complex[1:] ]
        cochain_complex += [ sparse.csc_matrix( (1, Bn.shape[1]), dtype=Bn.dtype) ] 

        # store the data members
        self.simplices = simplices
        self._chain_complex   = chain_complex
        self._cochain_complex = cochain_complex

    def complex_dimension(self):
        return len(self.cochain_complex()) - 1

    def complex(self):
        return self.simplices

    def chain_complex(self):
        return self._chain_complex

    def cochain_complex(self):
        return self._cochain_complex
        

