__all__ = ['regular_cube_complex']

from scipy import ones, ndim
import scipy

from pydec.mesh import regular_cube_mesh
from .cube_array import cube_array_boundary

class regular_cube_complex(list):
    """
    Represents the complex for a regular_cube_mesh
    """
    class complex_data:
        pass

    def __init__(self,mesh):
        if not isinstance(mesh,regular_cube_mesh):
            raise ValueError('expected a regular_cube_mesh')

        self.mesh = mesh

        self.__construct_hierarchy()
        
        self.vertices = self[0].cube_array.astype(float)

        #self.__test()

    def __repr__(self):
        output = ""
        output += "regular_cube_complex:\n"
        output += "  Shape:   " + str(self.mesh.bitmap.shape) + "\n"
        output += "  Complex:\n"
        for i in reversed(range(len(self))):
            output += "   %10d: %2d-D cubes\n" % (self[i].cube_array.shape[0],i)
        return output

    def __test(self):
        #test boundary operator
        for prev,next in zip(self[:-1],self[1:]):
            assert((prev.boundary*next.boundary).nnz == 0)

    def chain_complex(self):
        return [lvl.boundary for lvl in self]

    def cochain_complex(self):
        return [lvl.boundary.T.tocsr() for lvl in self[1:]] + \
                [scipy.sparse.csr_matrix((1,self[-1].cube_array.shape[0]),dtype='int8')]

    def complex_dimension(self):
        return self.mesh.dimension()
    
    def embedding_dimension(self):
        return self.mesh.dimension()

    def __construct_hierarchy(self):
        for i in range(self.complex_dimension() + 1):
            self.append(self.complex_data())

        self[-1].cube_array = self.mesh.cube_array()

        for i in reversed(range(self.complex_dimension())):
            faces,boundary = cube_array_boundary(self[i+1].cube_array,i+1)
            self[i  ].cube_array = faces
            self[i+1].boundary   = boundary

        self[0].boundary = scipy.sparse.csr_matrix((1,self[0].cube_array.shape[0]),dtype='int8')
