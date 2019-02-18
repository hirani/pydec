from numpy import vstack, hstack, ones, arange, empty, alltrue, array, \
        lexsort, ravel, hsplit, empty, ascontiguousarray, ndim

from scipy.sparse import csr_matrix

__all__ = ['cube_array_search', 'cube_array_boundary']

def cube_array_search(k_face_array,k_faces):
    """
    Find the row indices (of s) corresponding to the
    cubes stored in the rows of cube array v.
    It is assumed that the rows of s are sorted in
    lexicographical order.

    Example:

      k_face_array = array([[0,0,0],[0,0,1],[0,1,0],[1,0,1]])
      k_faces = array([[0,1,0],[0,0,1]])
      cube_array_searchsorted(k_face_array,k_faces)

    Returns:

      array([2,1])

    """

    if ndim(k_face_array) != 2 or ndim(k_faces) != 2:
        raise ValueError('expected rank 2 arrays')
        
    if k_face_array.shape[1] != k_faces.shape[1]:
        raise ValueError('number of columns must agree')

    # a dense array used to lookup k_face_array row indices 
    lookup_grid_dimensions = k_face_array.max(axis=0) + 1
    
    lookup_grid = empty(lookup_grid_dimensions,dtype=k_faces.dtype)
    lookup_grid[:] = -1
    lookup_grid[hsplit(k_face_array,k_face_array.shape[1])] = arange(k_face_array.shape[0],dtype=k_faces.dtype).reshape((-1,1))
    row_indices = lookup_grid[hsplit(k_faces,k_faces.shape[1])].reshape((-1))

    return row_indices


def cube_array_boundary(cubes,dim):
    """
    Compute the faces and boundary operator for a regular grid of N-cubes
    """
    cube_dim = dim
    top_dim  = cubes.shape[1] - cube_dim
    num_cubes = cubes.shape[0]
    num_faces = 2*num_cubes*cube_dim

    faces = empty((num_faces,cubes.shape[1]+1),dtype='int32')

    # Use simplex orientation to determine sign of boundary faces
    for i in range(cube_dim):
        # Faces originating at cube origin
        rows = faces[(2*i+0)*num_cubes:(2*i+1)*num_cubes]
        rows[:,:top_dim+i ]    = cubes[:,:top_dim+i   ]
        rows[:, top_dim+i:-2]  = cubes[:, top_dim+i+1:]
        rows[:,-2]             = arange(num_cubes)
        rows[:,-1]             = (-1)**(i+1)

        # Faces originating at other corners of cube
        rows = faces[(2*i+1)*num_cubes:(2*i+2)*num_cubes]
        rows[:,:top_dim+i ]    = cubes[:,:top_dim+i   ]
        rows[:, top_dim+i:-2]  = cubes[:, top_dim+i+1:]
        rows[:,-2]             = arange(num_cubes)
        rows[:,-1]             = (-1)**(i+2)
        rows[arange(rows.shape[0]),cubes[:,top_dim+i]] += 1

    #sort rows
    faces = faces[lexsort([faces[:,i] for i in reversed(range(faces.shape[1]-2))])]

    #find unique faces
    face_mask = -hstack((array([False]),alltrue(faces[1:,:-2] == faces[:-1,:-2],axis=1)))

    unique_faces = faces[face_mask,:-2].copy()

    #compute CSR representation for boundary operator
    indptr  = hstack((arange(num_faces)[face_mask],array([num_faces])))
    indices = ascontiguousarray(faces[:,-2])
    data    = ascontiguousarray(faces[:,-1].astype('int8'))

    shape = (len(unique_faces),num_cubes)
    boundary_operator = csr_matrix((data,indices,indptr), shape )

    return unique_faces,boundary_operator
