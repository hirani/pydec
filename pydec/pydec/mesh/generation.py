__all__ = ['simplicial_grid_2d','cube_grid']


from scipy import zeros,resize,arange,ravel,concatenate,matrix, \
     transpose,prod,mgrid,ndindex,sum,array,cumprod,tile,ones
import scipy


def simplicial_grid_2d(n):
    """
    Create an NxN 2d grid in the unit square
    
    The number of vertices along each axis is (N+1) for a total of (N+1)x(N+1) vertices
    
    A tuple (vertices,indices) of arrays is returned
    """
    vertices = zeros(((n+1)**2,2))
    vertices[:,0] = ravel(resize(arange(n+1),(n+1,n+1)))
    vertices[:,1] = ravel(transpose(resize(arange(n+1),(n+1,n+1))))
    vertices /= n
    
    indices = zeros((2*(n**2),3),scipy.int32)

    
    t1 = transpose(concatenate((matrix(arange(n)),matrix(arange(1,n+1)),matrix(arange(n+2,2*n+2))),axis=0))
    t2 = transpose(concatenate((matrix(arange(n)),matrix(arange(n+2,2*n+2)),matrix(arange(n+1,2*n+1))),axis=0))
    first_row = concatenate((t1,t2))
    
    for i in xrange(n):       
        indices[(2*n*i):(2*n*(i+1)),:] = first_row + i*(n+1)
    
    return (vertices,indices)


def cube_grid(dims):
    """
    Return a regular nD-cube mesh with given shape.

    Eg.
      cube_grid_nd((2,2))   -> 2x2   - 2d mesh (x,y)
      cube_grid_nd((4,3,2)) -> 4x3x2 - 3d mesh (x,y,z)

    Eg.
    
      v,i = cube_grid_nd((2,1))

      v =
      array([[ 0.,  0.],
             [ 1.,  0.],
             [ 2.,  0.],
             [ 0.,  1.],
             [ 1.,  1.],
             [ 2.,  1.]])

      i = 
      array([[[0, 3],
              [1, 4]],

             [[1, 4],
              [2, 5]]])

    """
    dims = tuple(dims)
    
    vert_dims = tuple(x+1 for x in dims)
    N = len(dims)
    
    vertices = zeros((prod(vert_dims),N))
    grid     = mgrid[tuple(slice(0,x,None) for x in reversed(vert_dims))]
    for i in range(N):
        vertices[:,i] = ravel(grid[N-i-1])


    #construct one cube to be tiled
    cube  = zeros((2,)*N,dtype='i')
    cycle = array([1] + list(cumprod(vert_dims)[:-1]),dtype='i')
    for i in ndindex(*((2,)*N)):
        cube[i] = sum(array(i) * cycle)
        cycle = array([1] + list(cumprod(vert_dims)[:-1]),dtype='i')


    #indices of all vertices which are the lower corner of a cube
    interior_indices = arange(prod(vert_dims)).reshape(tuple(reversed(vert_dims))).T
    interior_indices = interior_indices[tuple(slice(0,x,None) for x in dims)]

    indices = tile(cube,(prod(dims),) + (1,)*N) + interior_indices.reshape((prod(dims),) + (1,)*N)
    
    return (vertices,indices)
