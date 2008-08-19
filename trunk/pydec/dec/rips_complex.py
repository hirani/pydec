from numpy import ones, arange, array, asarray, hstack, empty, lexsort
from scipy.sparse import csr_matrix, csc_matrix

from pydec.mesh.simplex import simplex
from pydec.math import kd_tree
from pydec.util import flatten

from simplex_array import simplex_array_searchsorted
 
class rips_complex:
    def __init__(self, vertices, r):
        vertices = asarray(vertices)
        
        tree = kd_tree(vertices)

        L = [tree.in_sphere(v,r) for v in vertices]

        i = arange(len(vertices)).repeat( [len(x) for x in L] )
        j = array( flatten(L) )
        edges = hstack((i.reshape(-1,1),j.reshape(-1,1)))

        simplices = rips_simplices(edges, num_nodes=vertices.shape[0], \
                k=vertices.shape[1])

        chain_complex = rips_chain_complex(simplices)


        self.vertices = vertices
        self.simplices = simplices
        self._chain_complex = chain_complex

        Bn = chain_complex[-1]
        self._cochain_complex  = [ B.T for B in chain_complex[1:] ]
        self._cochain_complex += [ csr_matrix( (1, Bn.shape[1]), dtype=Bn.dtype) ] 
    
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
    Bs = []

    dtype = 'int8'

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



def rips_simplices(edges,num_nodes=None,k=2):
    edges = asarray(edges,dtype='int32')

    if num_nodes is None:
        #infer from edge array
        num_nodes = edges.max() + 1

    simplices = []
        
    if len(edges.shape) != 2 and edges.shape[1] != 2:
        raise ValueError,'expected N by 2 array for argument edges'

    # node simplices 
    simplices.append( arange( num_nodes, dtype=edges.dtype).reshape(-1,1) )

    if k == 0 or len(edges) == 0: return simplices
    
    edges = edges[edges[:,0] < edges[:,1]]   # orient edges 
    edges = edges[lexsort( edges.T[::-1] )]  # sort edges 

    simplices.append( edges )

    if k == 1: return simplices

    # E[i,j] == 1 iff (i,j) is an edge and i < j
    E = csr_matrix( (ones(len(edges),dtype='int8'), edges.T),shape=(num_nodes,num_nodes))

    k_faces = edges

    for n in range(1,k):
        indptr  = arange(0, (n+1) * (len(k_faces)+1), (n+1) )
        indices = k_faces.ravel()
        data    = ones( len(indices), dtype='int8')
          
        F = csr_matrix((data,indices,indptr),shape=(len(k_faces),num_nodes))

        FE = F*E
        FE.data[FE.data != n+1] = 0
        FE.eliminate_zeros()
        FE.sort_indices()

        row = FE.tocoo(copy=False).row
        k_faces = hstack((k_faces[row],FE.indices.reshape(-1,1)))

        if len(k_faces) == 0:
            return simplices

        simplices.append( k_faces )

    return simplices



