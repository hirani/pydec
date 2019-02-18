__all__ = ['Simplex','SimplicialMesh','simplex','simplicial_mesh']

from pydec.math import signed_volume,relative_parity,combinations
from .base_mesh import base_mesh


from scipy import asarray
import numpy,scipy


class simplex(tuple):
    def __new__(cls,s,parity=0):
        obj = tuple.__new__(cls, sorted(s))
        obj.parity =  relative_parity(obj,s) ^ parity
        return obj
    
    def __repr__(self):
        return 'simplex(' + tuple.__repr__(self) + ',parity=' + str(self.parity) + ')'
    
    def boundary(self):
        """
        A list of oriented simplicies in the boundary of this simplex
        """
        return [ simplex(self[:n] + self[n+1:], (self.parity + n) % 2) for n in range(len(self)) ]



class simplicial_mesh(base_mesh):
    """Simplicial mesh

    Can be instantiated in several ways:
        - simplicial_mesh(V,E)
            - where V and E are arrays of vertices and simplex indices
        - simplicial_mesh( D )
            - where D is a dictionary with keys 'vertices' and 'elements'
    

    Examples
    ========

    >>> from numpy import array
    >>> from pydec.mesh import simplicial_mesh
    >>> V = array([[0,0],[1,0],[0,1]]) # mesh vertices
    >>> E = array([[0,1,2]])           # mesh indices
    >>> simplicial_mesh(V,E)

    >>> D = {'vertices' : V, 'elements' : E}
    >>> simplicial_mesh(data) 
    
    """

    def __init__(self,*args,**kwargs):

        if len(args) == 2:
            #try to parse as (vertices,indices)
            V,I = args
            self['vertices'] = asarray(V)
            self['elements'] = asarray(I)
        elif len(kwargs) == 2:
            self['vertices'] = asarray(kwargs['vertices'])
            self['elements'] = asarray(kwargs['indices'])
        elif len(args) == 1 and isinstance(args[0],dict):
            base_mesh.update(self,args[0])
        else:
            raise ValueError('unrecognized arguments')

        if numpy.ndim(self['elements']) != 2 or numpy.ndim(self['vertices']) != 2:
            raise ValueError('index and vertex arrays must have rank 2')

        if self['elements'].min() < 0 or self['elements'].max() > self['vertices'].shape[0]:
            raise ValueError('invalid index value')
            

    def __getattr__(self, attr):
        if attr == 'vertices':
            return self['vertices']
        elif attr in ['indices', 'elements']:        
            return self['elements']

        return base_mesh.__getattr__(self, attr)

    def __setattr__(self,attr,value):
        if attr == 'vertices':
            self['vertices'] = value
        elif attr in ['indices','elements']:        
            self['elements'] = value
        else:
            return base_mesh.__setattr__(self,attr,value)

    def __repr__(self):
        output = ""
        output += "simplicial_mesh< " + str(self.manifold_dimension()) + "D manifold, "
        output += str(self.embedding_dimension()) + "D embedding, "
        output += str(self['vertices'].shape[0]) + " vertices, "
        output += str(self['elements'].shape[0])  + " elements >\n"

        format_str = '\t%-16s %16s %10s\n'

        output +=  format_str % ('Data Names'.center(16),'Shape'.center(16),'Size (KB)'.center(16))
        for k,v in self.iteritems():
            output += format_str % (k,str(v.shape),str(v.nbytes/1024))
        return output      
    
    def manifold_dimension(self):
        return self['elements'].shape[1] - 1
    
    def embedding_dimension(self):
        return self['vertices'].shape[1]
    
    def boundary(self):        
        """
        Return a set() of the boundary simplices, i.e. the faces 
        of the top level simplices that occur only once
        """
        boundary_set = set()
        
        for row in self['elements']:
            s = simplex(row)     
            for b in s.boundary():
                if b in boundary_set:
                    boundary_set.remove(b) #b has occured twice
                else:
                    boundary_set.add(b)    #first occurance of b
        
        return boundary_set
    
    def skeleton(self,p):
        """
        Returns the p-skeleton (all the p-faces) of the mesh as a set()
        """        
        assert(0 <= p <= self.manifold_dimension())
        
        skeleton_set = set()
        
        for row in self.indices:
            for b in combinations(row,p+1):
                skeleton_set.add(simplex(b))
        
        return skeleton_set
    

    def orient(self):
        """
        Orient this SimplicialMesh. If the  manifold is of the same dimension as the
        embedding (e.g. triangle mesh in 2D, tet mesh in 3D) then the resultant mesh
        will be oriented so that the simplices have positive volume.
        
        If the mesh is not orientable an Exception is raised.
        """

        
        if self.manifold_dimension() == 0:
            #0-dimensional manifold is alway oriented
            return
        
        if self.manifold_dimension() == self.embedding_dimension():
            #orient w/ positive volumes
            num_flips = 0
            elements = self['elements']
            vertices = self['vertices']
            
            for row in elements:
                pts = vertices[row]
                if signed_volume(pts) < 0:
                    num_flips += 1
                    temp = row[0]
                    row[0] = row[1]
                    row[1] = temp
            print("Flipped", num_flips,"simplices")
            return

        raise NotImplementedError
    
        simplex_to_index = {}
        for index,row in enumerate(self['elements']):
            simplex_to_index[simplex(row)] = index
        
        faces = self.skeleton(self.manifold_dimension() - 1)
        face_to_simplex = dict.fromkeys(faces,set())
        for simplex,index in simplex_to_index.iterkeys():
            for b in simplex.boundary():
                face_to_simplex[b].add(index)
            
        simplex_neigbors = [[]]*len(self['elements'])
        for simplex,index in simplex_to_index.iterkeys():            
            for b in simplex.boundary():
                simplex_neigbors[index].append(face_to_simplex[b] - set([index]))
        
        print(simplex_neighbors)
                

#for backwards compatibility
Simplex = simplex
SimplicialMesh = simplicial_mesh
