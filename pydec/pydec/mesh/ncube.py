__all__ = ['nCube','nCubeMesh','RegularCubeMesh']


from numpy import ndarray,array,asarray,ndim,bitwise_xor,eye,hstack,vstack,arange,zeros
from .simplex import Simplex


class RegularCubeMesh:
    """
    A regular grid of hypercubes.


    Examples:
    
       # create a 2x2 cube mesh
       bitmap = ones((2,2),dtype='bool')
       c_mesh = RegularCubeMesh(bitmap)

       # creates a 3x3 cube mesh with a center hole
       bitmap = ones((3,3),dtype='bool')
       bitmap[1,1] = False
       c_mesh = RegularCubeMesh(bitmap)

       # creates a 10x10x10 cube mesh with a center hole
       bitmap = ones((10,10,10),dtype='bool')
       bitmap[5,5,5] = False
       c_mesh = RegularCubeMesh(bitmap)
    """
    def __init__(self,bitmap):
        self.bitmap = asarray(bitmap,dtype='bool')

    def cube_array(self):
        """
        Return a cube array that represents this mesh's bitmap
        """
        cubes = vstack(self.bitmap.nonzero()).transpose()
        applied_zeroes = zeros((cubes.shape[0], ndim(self.bitmap)), dtype=cubes.dtype)
        cubes = hstack((cubes, applied_zeroes + arange(ndim(self.bitmap))))

        return cubes
    
    def dimension(self):
        return ndim(self.bitmap)




class nCube:
    def __init__(self,s,dtype='int32'):
        self.indices = array(s,dtype=dtype)
        assert(self.indices.shape == (1,) or self.indices.shape == (2,)*ndim(self.indices))
        
        self.compute_corner_simplex()

    def compute_corner_simplex(self):
        if ndim(self.indices) < 2:
            self.corner_simplex = Simplex(self.indices)
        else:
            corner_value = self.indices.min()
            corner_index = (self.indices == corner_value).nonzero()

            rest = self.indices[[tuple(x) for x in bitwise_xor(eye(ndim(self.indices),dtype=int),array(corner_index))]]

            parity = sum(corner_index)[0] % 2

            self.corner_simplex = Simplex([corner_value] + rest.tolist(),parity)
        
    def __str__(self):
        return 'nCube(' + ndarray.__str__(self.indices) + ')'
    
    def __hash__(self):
        return hash(self.corner_simplex)
    
    def __eq__(self,other):
        return self.corner_simplex == other

## Ideas for boundary()
##In [17]: A.take([0],axis=0)
##Out[17]: 
##array([[[0, 1],
##        [2, 3]]])
##
##In [18]: A.take([1],axis=0)
##Out[18]: 
##array([[[4, 5],
##        [6, 7]]])
    def boundary(self):
        raise NotImplementedError


    def relative_orientation(self,other):
        """
        Determine whether two cubes that represent the same
        face have the same orientation or opposite orientations

          Returns:
             False if same orientation
             True if opposite orientation
        """
        if self.corner_simplex != other.corner_simplex:
            raise ValueError('Cubes do not share the same vertices')
        return self.corner_simplex.parity ^ other.corner_simplex.parity




        
    
class nCubeMesh:
    def __init__(self,indices,vertices):
        self.indices  = asarray(simplices,dtype='i')
        self.vertices = asarray(vertices, dtype='d')

    def manifold_dimension(self):
        if ndim(self.indices) >= 2:
            return ndim(self.indices)
        else:
            return self.indices.shape[1]
    
    def embedding_dimension(self):
        return self.vertices.shape[1]



