__all__ = ['regular_cube_mesh']

from numpy import array,asarray,hstack,vstack,arange,zeros,ndim


class regular_cube_mesh:
    """
    A regular grid of hypercubes.


    Examples:
    
       # create a 2x2 cube mesh
       bitmap = ones((2,2),dtype='bool')
       c_mesh = regular_cube_mesh(bitmap)

       # creates a 3x3 cube mesh with a center hole
       bitmap = ones((3,3),dtype='bool')
       bitmap[1,1] = False
       c_mesh = regular_cube_mesh(bitmap)

       # creates a 10x10x10 cube mesh with a center hole
       bitmap = ones((10,10,10),dtype='bool')
       bitmap[5,5,5] = False
       c_mesh = regular_cube_mesh(bitmap)
    """
    def __init__(self,bitmap):
        self.bitmap = asarray(bitmap,dtype='bool')

    def cube_array(self):
        """
        Return a cube array that represents this mesh's bitmap
        """
        cubes = vstack(self.bitmap.nonzero()).transpose().astype('int32')
        cubes = hstack((cubes,zeros((cubes.shape[0],ndim(self.bitmap)),dtype=cubes.dtype) + arange(ndim(self.bitmap),dtype=cubes.dtype)))

        return cubes
    
    def dimension(self):
        return ndim(self.bitmap)
