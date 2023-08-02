__all__ = ['simplex_array_searchsorted','simplex_array_boundary','simplex_array_parity']


from scipy import ravel, zeros, ones, arange, empty, alltrue, array, lexsort, \
    hstack, vstack, ndim, bincount, cumsum, ascontiguousarray, zeros_like, \
    concatenate, asarray
from scipy.sparse import csr_matrix

def simplex_array_searchsorted(s, v):
    """Find the row indices (of s) corresponding to the simplices stored 
    in the rows of simplex array v.  The rows of s must be stored in 
    lexicographical order.

    Example
    -------

    >>> from numpy import array
    >>> s = array([[0,1],[0,2],[1,2],[1,3]])
    >>> v = array([[1,2],[0,2]])
    >>> simplex_array_searchsorted(s,v)
    array([2, 1])

    """

    s = asarray(s)
    v = asarray(v)

    if ndim(s) != 2 or ndim(v) != 2:
        raise ValueError('expected rank 2 arrays')

    if s.shape[1] != v.shape[1]:
        raise ValueError('number of columns must agree')
   
    # compute row indices by sorting both arrays together
    Ns = s.shape[0]
    Nv = v.shape[0]
    
    perm = lexsort(vstack((s,v))[:,::-1].T)
    
    flags = concatenate( (ones(Ns,dtype=int),zeros(Nv,dtype=int)) )
    indices = empty(Ns+Nv, dtype=int)
    indices[perm] = cumsum(flags[perm])
    indices = indices[Ns:].copy()
    indices -= 1

    return indices




def simplex_array_parity(s):
    """Compute the relative parity of an array of simplices
    """
    s = s.copy()

    M,N = s.shape

    # number of transpositions used to sort the
    # indices of each simplex (row of s)
    trans = zeros_like(s[:,0])
    seq   = arange(M)

    for i in range(N-1):
        pos = s.argmin(axis=1)
        s[seq,pos] = s[:,0]
        pos.clip(0,1,pos)
        trans += pos
        s = s[:,1:]

    trans %= 2 #compute parity

    return trans




def simplex_array_boundary(s,parity):
    """
    Compute the boundary faces and boundary operator of an
    array of simplices with given simplex parities

    E.g.
    
      For a mesh with two triangles [0,1,2] and [1,3,2], the second
      triangle has opposite parity relative to sorted order.
      
      simplex_array_boundary(array([[0,1,2],[1,2,3]]),array([0,1]))
      
    """
    #TODO handle edge case as special case
    
    num_simplices     = s.shape[0]
    faces_per_simplex = s.shape[1]
    num_faces         = num_simplices * faces_per_simplex

    orientations = 1 - 2*parity

    #faces[:,:-2] are the indices of the faces
    #faces[:,-2]  is the index of the simplex whose boundary produced the face
    #faces[:,-1]  is the orientation of the face in the boundary of the simplex
    faces = empty((num_faces,s.shape[1]+1),dtype=s.dtype)
    for i in range(faces_per_simplex):
        rows = faces[num_simplices*i:num_simplices*(i+1)]

        rows[:,  : i] = s[:,   :i]
        rows[:,i :-2] = s[:,i+1: ]
        rows[:, -2  ] = arange(num_simplices)
        rows[:, -1  ] = ((-1)**i)*orientations

    #sort rows
    faces = faces[lexsort( faces[:,:-2].T[::-1] )]

    #find unique faces
    face_mask    = ~hstack((array([False]), alltrue(faces[1:,:-2] == faces[:-1,:-2], axis=1)))

    unique_faces = faces[face_mask,:-2]

    #compute CSR representation for boundary operator
    csr_indptr  = hstack((arange(num_faces)[face_mask],array([num_faces])))
    csr_indices = ascontiguousarray(faces[:,-2])
    csr_data    = faces[:,-1].astype('int8')
  
    shape = (len(unique_faces),num_simplices)   
    boundary_operator = csr_matrix((csr_data,csr_indices,csr_indptr), shape)

    return unique_faces,boundary_operator





####################################    
## Fast C implementations
####################################
#
#
#import scipy
#
#def simplex_array_searchsorted(s,v):
#    """
#    Find the row indices (of s) corresponding to the
#    simplices stored in the rows of simplex array v.
#    It is assumed that the rows of s are sorted in
#    lexicographical order.
#
#    Example:
#
#      s = array([[0,1],[0,2],[1,2],[1,3]])
#      v = array([[1,2],[0,2]])
#      simplex_array_searchsorted(s,v)
#
#    Returns:
#
#      array([2,1])
#
#    """
#
#    if ndim(s) != 2 or ndim(v) != 2:
#        raise ValueError,'expected rank 2 arrays'
#
#    if s.shape[1] != v.shape[1]:
#        raise ValueError,'number of columns must agree'
#
#    s_row = s.shape[0]
#    v_row = v.shape[0]
#    s_col = s.shape[1]
#    s_max = int(s[:,0].max())
#    
#    first_index = cumsum(hstack((array([0]),bincount(s[:,0]))))
#    indices = empty(v.shape[0],dtype=s.dtype)
#
#    code = """
#    #line 45 "simplex_array.py"
#
#    for(int i = 0; i < v_row; i++){
#
#      int v_i0 = v(i,0);
#        
#      if (v(i,0) > s_max)
#        py::fail(PyExc_ValueError, "index exceeds expected range");
#    
#      int row_start = first_index(v_i0);
#      int row_end   = first_index(v_i0+1);
#
#      int row_index = -1;
#      for(int k = row_start; k < row_end; k++){
#        bool mismatch = false;
#        for(int j = 1; j < s_col; j++){
#          if(v(i,j) != s(k,j)){
#            mismatch = true;
#            break;
#          }
#        }
#        if(!mismatch){
#          row_index = k;
#          break;
#        }
#      }
#      if (row_index == -1)
#        py::fail(PyExc_ValueError, "simplex not found");
#
#      indices(i) = row_index;  
#    }
#    """
#
#    err = scipy.weave.inline(code,
#                             ['s','v','first_index', 'indices', 's_row', 'v_row', 's_col','s_max'],
#                             type_converters = scipy.weave.converters.blitz,
#                             compiler = 'gcc')
#
#    return indices
#
#
#
#def simplex_array_parity(s):
#    """
#    Compute the relative parity of an array of simplices
#    """
#    perms = s.argsort()
#
#
#    n_rows,n_cols = perms.shape
#
#    parity = zeros(n_rows,dtype=perms.dtype)
#    seen   = zeros(n_cols,dtype=perms.dtype)
#    
#    code = """
#        #line 26 "relaxation.py"
#
#        for(int i = 0; i < n_rows; i++){
#           int num_cycles = 0;
#           for(int j = 0; j < n_cols; j++){
#             
#             if(seen(j) == 1) continue;
#
#             int k = j;
#             
#             num_cycles++;
#             while ( true ){
#               seen(k) = 1;
#               k = perms(i,k);
#               if (j == k) break;
#             }
#           }
#           for(int j = 0; j < n_cols; j++){ seen(j) = 0; } //reset seen
#           parity(i) = (n_cols - num_cycles) % 2;
#        }
#        """
#
#    from time import clock
#
#
#    # compiler keyword only needed on windows with MSVC installed
#    err = scipy.weave.inline(code,
#                       ['perms', 'parity', 'seen', 'n_rows', 'n_cols' ],
#                       type_converters = scipy.weave.converters.blitz,
#                       compiler = 'gcc')
#
#    return parity


