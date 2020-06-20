__all__ = ['SimplicialComplex','simplicial_complex']

from warnings import warn

import numpy
import scipy
from numpy import array, dot, inner, ones, cross, copysign
from scipy import sparse, zeros, asarray, mat, hstack

import pydec
from pydec.mesh.simplex import simplex, simplicial_mesh
from pydec.math import (circumcenter, unsigned_volume, signed_volume,
                            circumcenter_barycentric)

from .simplex_array import simplex_array_parity, simplex_array_boundary
from .cochain import cochain


class simplicial_complex(list):
    """simplicial complex

    This can be instantiated in several ways:
        - simplicial_complex( (V,S) )
            - where V and S are arrays of vertices and simplex indices
        - simplicial_complex( M )
            - where M is a simplicial_mesh object


    Examples
    ========
    >>> from pydec import simplicial_complex, simplicial_mesh
    >>> V = [[0,0],[1,0],[1,1],[0,1]]
    >>> S = [[0,1,3],[1,2,3]]
    >>> M = simplicial_mesh(V, S)
    >>> simplicial_complex( (V,S) )
    >>> simplicial_complex( M )


    """

        
    def __init__(self, arg1, arg2=None):                    

        if arg2 is not None:
            warn('initializing a simplicial_complex with' \
                    ' format (vertices,simplices) is deprecated', \
                    DeprecationWarning)
            arg1 = (arg1,arg2)

        if isinstance(arg1,simplicial_mesh):
            self.mesh = arg1
        elif isinstance(arg1,tuple):
            self.mesh = simplicial_mesh(arg1[0], arg1[1])
        else:
            raise ValueError('unrecognized constructor usage') 

        #TODO is this necessary?
        #if self.mesh['vertices'].dtype != 'float32':
        #    self.mesh['vertices'] = asarray(self.mesh['vertices'],dtype='float64')

        self.vertices  = self.mesh['vertices']
        self.simplices = self.mesh['elements']
    
        self.build_complex(self.simplices)
    
    def __repr__(self):
        output = ""
        output += "simplicial_complex:\n"
        output += "  complex:\n"
        for i in reversed(range(len(self))):
            output += "   %10d: %2d-simplices\n" % (self[i].num_simplices,i)
        return output

        
    def complex_dimension(self):
        return self.simplices.shape[1] - 1

    def embedding_dimension(self):
        return self.vertices.shape[1]

    def chain_complex(self):
        return [x.boundary for x in self]
        
    def cochain_complex(self):
        return [x.d for x in self]

    def complex(self):
        return [x.simplices for x in self]

    def get_cochain(self, dimension, is_primal=True):
        N = self.complex_dimension()

        if not 0 <= dimension <= N:
            raise ValueError('invalid dimension (%d)' % dimension)
        
        c   = cochain(self, dimension, is_primal)        
        
        if is_primal:
            c.v = zeros(self[dimension].num_simplices)
        else:
            c.v = zeros(self[N - dimension].num_simplices)
        
        return c   
    
    def get_cochain_basis(self, dimension, is_primal=True):
        N = self.complex_dimension()

        if not 0 <= dimension <= N:
            raise ValueError('invalid dimension (%d)' % dimension)
        
        c   = cochain(self, dimension, is_primal)        
        
        if is_primal:
            c.v = sparse.identity(self[dimension].num_simplices)
        else:
            c.v = sparse.identity(self[N - dimension].num_simplices)
        
        return c   

    def build_complex(self, simplex_array):
        """Compute faces and boundary operators for all dimensions"""
        N,K = simplex_array.shape

        s = simplex_array.copy()
        parity = simplex_array_parity(s)
        s.sort()

        simplices = [s]
        chain_complex = []
        parities = [parity]

        while s.shape[1] > 1:
            s,boundary = simplex_array_boundary(s,parity)
            parity = zeros(s.shape[0],dtype=s.dtype)
            
            simplices.append( s )
            chain_complex.append( boundary )
            parities.append( parity )

        B0 = sparse.csr_matrix( (1,len(s)), dtype='uint8')
        chain_complex.append( B0 )

        simplices     = simplices[::-1]
        chain_complex = chain_complex[::-1]
        parities      = parities[::-1]

        Bn = chain_complex[-1]
        cochain_complex  = [ B.T for B in chain_complex[1:] ]
        cochain_complex += [ sparse.csc_matrix( (1, Bn.shape[1]), dtype=Bn.dtype) ] 

        for n in range(K):
            data = self.data_cache()

            data.d              = cochain_complex[n]
            data.boundary       = chain_complex[n]
            data.complex        = self
            data.dim            = n
            data.simplices      = simplices[n]
            data.num_simplices  = len(data.simplices)
            data.simplex_parity = parities[n]

            self.append(data)
     
        

    def construct_hodge(self):
        """Construct the covolume Hodge star for all levels"""        
        
        for dim,data in enumerate(self):
            form_size = data.num_simplices
            data.star = sparse.lil_matrix((form_size,form_size))                      
            data.star_inv = sparse.lil_matrix((form_size,form_size))
            
            stardiag = data.dual_volume / data.primal_volume
            N = len(stardiag)
            data.star     = sparse.spdiags([stardiag],    [0], N, N, format='csr')
            data.star_inv = sparse.spdiags([1.0/stardiag],[0], N, N, format='csr')
            
            #Choose sign of star_inv to satisfy (star_inv * star) = -1 ^(k*(n-k))
            n,k = self.complex_dimension(),dim            
            data.star_inv *= (-1) ** (k * (n - k))                   

    def compute_bary_circumcenters(self,dim):
        """Compute circumcenters for all simplices at a given
        dimension in barycentric coordinates
        """
        data = self[dim]            
        data.bary_circumcenter = zeros((data.num_simplices,
                                        dim + 1))
        
        for i,s in enumerate(data.simplices):
            pts = self.vertices[[x for x in s],:]
            data.bary_circumcenter[i] = circumcenter_barycentric(pts)
            
    def compute_circumcenters(self,dim):
        """Compute circumcenters for all simplices at a given dimension
        """
        data = self[dim]            
        data.circumcenter = zeros((data.num_simplices,
                                       self.embedding_dimension()))
        
        for i,s in enumerate(data.simplices):
            pts = self.vertices[[x for x in s],:]
            data.circumcenter[i] = circumcenter(pts)[0]

    def compute_primal_volume(self,dim):
        """Compute the volume of all simplices for a given dimension

        If the top simplex is of the same dimension as its embedding,
        the signed volume is computed.        
        """
        data = self[dim]
        data.primal_volume = zeros((data.num_simplices,))

        if dim == self.embedding_dimension():
            for i,s in enumerate(self.simplices):
                pts = self.vertices[s,:]            
                data.primal_volume[i] = signed_volume(pts)
        else:

            for i,s in enumerate(data.simplices):
                pts = self.vertices[s,:]            
                data.primal_volume[i] = unsigned_volume(pts)

    def compute_dual_volume(self):
        """Compute dual volumes for simplices of all dimensions
        """
        for dim,data in enumerate(self):
            data.dual_volume = zeros((data.num_simplices,))

        temp_centers = zeros((self.complex_dimension() + 1,
                                  self.embedding_dimension()))
        temp_bcenters = [zeros(dim) for dim in
                             range(1, self.complex_dimension() + 2)]
        temp_signs = ones(self.complex_dimension() + 1)
        for i,s in enumerate(self.simplices):
            self.__compute_dual_volume(simplex(s),
                                       None,
                                       temp_centers,
                                       temp_bcenters,
                                       temp_signs,
                                       self.complex_dimension())
            
    def __compute_dual_volume(self, s, parent, pts, bpts, signs, dim):    
        ## s: simplex whose dual volume is desired (NumPy array)
        ## parent: simplex initiating this call whose face is s (NumPy
        ## array) 
        ## pts: circumcenters at all dimensions (NumPy array)
        ## bpts: circumcenters at all dimensions in barycentric 
        ## coordinates (Python list)

        ## Remark: Sign s is determined by location of cicumcenter
        ## with respect to a face of a simplex and the opposite vertex
        ## of the face in the simplex. It is +1 if the circumcenter
        ## and the opposite vertex lie on the same side with respect
        ## to the halfplane generated by the face. It is -1 if they
        ## lie on opposite sides of this halfplane. Consequently, sign
        ## can be determined as the sign of that component of the
        ## barycentric coordiante of the circumcenter which
        ## corresponds to the opposite vertex.

        data = self[dim]
        index = data.simplex_to_index[s]  
        pts[dim] = data.circumcenter[index]
        bpts[dim] = data.bary_circumcenter[index]
        sgn = 1
        if parent is not None:
            opposite_vertex = list(set(parent) - set(s))[0]
            ov_index = list(parent).index(opposite_vertex)
            signs[dim] = copysign(1, bpts[-1][ov_index])
        for i in range(dim, self.complex_dimension()):
            sgn *= signs[i]
        data.dual_volume[index] += sgn * unsigned_volume(pts[dim:,:])
        if dim > 0:
            for bs in s.boundary():
                self.__compute_dual_volume(bs, s, pts, bpts, signs, dim-1)


    def boundary(self):
        """Return a list of the boundary simplices, i.e. the faces of the top level simplices that occur only once
        """
    
        assert(self.complex_dimension() > 0) 
        old_keys = self[self.complex_dimension() - 1].simplex_to_index.keys()
        face_count = dict.fromkeys(old_keys, 0)
        
        for s in self[self.complex_dimension()].simplex_to_index.keys():
            for f in s.boundary():
                face_count[f] += 1
        
        boundary_faces = [f for f, count in face_count.items() if count == 1]
        return boundary_faces
        

    class data_cache:
        """caches the result of costly operations"""
        def __getattr__(self,attr):
            #print "constructing: ",attr
            if attr == "d":
                return self.d.tocsr()
            elif attr == "star":
                self.complex.construct_hodge()
                return self.star
            elif attr == "star_inv":
                self.complex.construct_hodge()
                return self.star_inv
            elif attr == "circumcenter":
                self.complex.compute_circumcenters(self.dim)
                return self.circumcenter
            elif attr == "bary_circumcenter":
                self.complex.compute_bary_circumcenters(self.dim)
                return self.bary_circumcenter
            elif attr == "primal_volume":
                self.complex.compute_primal_volume(self.dim)
                return self.primal_volume
            elif attr == "dual_volume":
                self.complex.compute_dual_volume()
                return self.dual_volume
            elif attr == "simplex_to_index":
                self.simplex_to_index = dict((simplex(x), i) for i, x in enumerate(self.simplices))
                return self.simplex_to_index
            elif attr == "index_to_simplex":
                self.simplex_to_index = dict((i, simplex(x)) for i, x in enumerate(self.simplices))
                return self.simplex_to_index
            else:
                raise AttributeError(attr + " not found")
    

#for backwards compatibility
SimplicialComplex = simplicial_complex
