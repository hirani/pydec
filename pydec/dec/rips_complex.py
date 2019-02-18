from numpy import ones, arange, array, asarray, hstack, empty, lexsort
from scipy.sparse import csr_matrix, csc_matrix

from pydec.mesh.simplex import simplex
from pydec.math import kd_tree
from pydec.util import flatten

from .simplex_array import simplex_array_searchsorted

__all__ = ['rips_complex']

class rips_complex:
    def __init__(self, vertices, delta):
        """Construct a Rips complex

        Construct a Rips for an array of vertices and a radius delta.

        Parameters
        ----------
        vertices : array_like
            An N-by-D array of vertex coordinates where N is the 
            number of vertices and D is the dimension of the space
        delta : float
            Radius used to determine the connectivity of the Rips complex.
            Vertices which are separated by no more than distance delta
            are connected by an edge in the 1-skeleton of the Rips complex.

        Examples
        --------
        >>> from pydec.dec import rips_complex
        >>> from numpy import array
        >>> vertices = array([[0,0],[2,0],[2,2],[0,2],[1,1]], dtype='float')
        >>> delta = 2.1
        >>> rc = rips_complex(vertices, delta)
        >>> print rc.simplices[0]
        [[0]
         [1]
         [2]
         [3]
         [4]]
        >>> print rc.simplices[1]
        [[0 1]
         [0 3]
         [0 4]
         [1 2]
         [1 4]
         [2 3]
         [2 4]
         [3 4]]
        >>> print rc.simplices[2]
        [[0 1 4]
         [0 3 4]
         [1 2 4]
         [2 3 4]]

        """

        vertices = asarray(vertices)
        
        # construct kD-Tree for spatial queries
        tree = kd_tree(vertices)

        # for each vertex, compute all vertices within distance r
        L = [tree.in_sphere(v, delta) for v in vertices]

        # convert list of list structure into array of edges
        i = arange(len(vertices)).repeat( [len(x) for x in L] )
        j = array( flatten(L) )
        edges = hstack((i.reshape(-1,1),j.reshape(-1,1)))

        # compute simplices of the Rips complex from edges (1-skeleton of Rips complex)
        simplices = rips_simplices(vertices.shape[0], edges, vertices.shape[1])

        # compute boundary operators of the Rips complex
        chain_complex = rips_chain_complex(simplices)

        # store the cochain complex
        Bn = chain_complex[-1]
        cochain_complex  = [ B.T for B in chain_complex[1:] ]
        cochain_complex += [ csr_matrix( (1, Bn.shape[1]), dtype=Bn.dtype) ] 
        
        # store the data members
        self.vertices  = vertices
        self.simplices = simplices
        self._chain_complex   = chain_complex
        self._cochain_complex = cochain_complex
    
    def complex_dimension(self):
        return len(self.cochain_complex()) - 1

    def embedding_dimension(self):
        return self.vertices.shape[1]

    def complex(self):
        return self.simplices

    def chain_complex(self):
        return self._chain_complex

    def cochain_complex(self):
        return self._cochain_complex
        

def rips_chain_complex(simplices):
    """Construct the boundary operators for the Rips complex 

    Constructs the boundary operators for the Rips complex for
    a given a list of simplex arrays.

    Parameters
    ----------
    simplices : list of arrays
        List of length D, where simplices[i] is the simplex
        array representing the i-dimensional simplices.

    Returns
    -------
    Bs : list of sparse matrices
        List of length D + 1 where Bs[i] is the boundary operator
        for the i-dimensional simplices.

    Examples
    --------
    TODO

    """
    Bs = []

    dtype = 'float32'

    simplices = [asarray(x) for x in simplices]

    B0 = csc_matrix( (1, len(simplices[0]) ), dtype=dtype )
    Bs.append( B0 )
    
    if len(simplices) == 1:
        B1 = csc_matrix( (len(simplices[0]), 1), dtype=dtype )
        Bs.append(B1)
        return Bs

    B1 = empty( 2*len(simplices[1]), dtype=dtype )
    B1[0::2] = -1
    B1[1::2] =  1
    B1 = csc_matrix( ( B1, simplices[1].ravel(), arange(0, 2*(len(simplices[1]) + 1), 2)), \
            shape=(len(simplices[0]), len(simplices[1])) )
    Bs.append( B1 )

    for f,s in zip(simplices[1:-1],simplices[2:]):
        k = f.shape[1]

        faces = empty(( len(s), (k+1), k),dtype=s.dtype)
        for i in range(k+1):
            faces[:,i,:i] = s[:,   : i ]
            faces[:,i,i:] = s[:,i+1:  ]

        indptr  = arange(0, (k+1)*(len(s) + 1), k+1)
        indices = simplex_array_searchsorted(f, faces.reshape(-1,k) )
        data    = empty( faces.shape[:-1], dtype=dtype )
        for i in range(k+1):
            data[:,i] = (-1)**i
        data = data.ravel()

        Bs.append( csc_matrix( (data,indices,indptr), shape=(len(f),len(s)) ) )

    return Bs



def rips_simplices(num_vertices, edges, k):
    """Compute simplices of the Rips complex

    Returns a list of simplex arrays representing
    the 0 to k-dimensional simplices of the Rips complex.

    Parameters
    ----------
    num_verticles : integer
        Number of vertices/points/nodes in the Rips complex
    edges : array_like
        An N-by-2 array containing the edges of the Rips complex
    k : integer
        Maximum simplex dimension

    Returns
    -------
    simplices : list of arrays
        List of length k+1, where simplices[i] is the simplex
        array representing the i-dimensional simplices.

    Examples
    --------
    >>> from numpy import array
    >>> edges = array([[0,1],[0,2],[1,0],[1,2],[2,0],[2,1]])
    >>> simplices = rips_simplices(3, edges, 2)
    >>> print simplices[0]
    [[0]
     [1]
     [2]]
    >>> print simplices[1]
    [[0 1]
     [0 2]
     [1 2]]
    >>> print simplices[2]
    [[0 1 2]]

    """

    edges = asarray(edges, dtype='int32')

    simplices = []
        
    if len(edges.shape) != 2 and edges.shape[1] != 2:
        raise ValueError('expected N by 2 array for argument edges')

    # node simplices 
    simplices.append( arange(num_vertices, dtype=edges.dtype).reshape(-1,1) )

    # no work to do
    if k == 0: return simplices
    
    edges = edges[edges[:,0] < edges[:,1]]   # orient edges 
    if edges.shape[0] != 0:
        edges = edges[lexsort( edges.T[::-1] )]  # sort edges 
        simplices.append( edges )

    if k == 1 or edges.shape[0] == 0: return simplices

    # E[i,j] == 1 iff (i,j) is an edge and i < j
    E = csr_matrix( (ones(len(edges), dtype='int8'), edges.T), shape=(num_vertices,num_vertices))

    k_faces = edges

    for n in range(1,k):
        if k_faces.shape[0] == 0:
            break

        indptr  = arange(0, (n+1) * (len(k_faces)+1), (n+1) )
        indices = k_faces.ravel()
        data    = ones(len(indices), dtype='int8')
          
        F = csr_matrix((data, indices, indptr), shape=(len(k_faces), num_vertices))

        # FE[i,j] == n+1 if all vertices in face i are adjacent to vertex j
        FE = F*E
        FE.data[FE.data != n+1] = 0
        FE.eliminate_zeros()
        FE.sort_indices()

        row = FE.tocoo(copy=False).row
        k_faces = hstack( (k_faces[row],FE.indices.reshape(-1,1)) )

        if len(k_faces) == 0:
            return simplices

        simplices.append( k_faces )

    return simplices

