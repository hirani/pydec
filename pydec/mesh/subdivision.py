__all__ = ['loop_subdivision','triangulate_ncube']

from pydec.math.combinatorial import combinations
from scipy import concatenate,matrix,ravel,array,ndim,vstack,hstack,zeros,arange,tile


def loop_subdivision(vertices,simplices):
    """
    Given a triangle mesh represented by the matrices (vertices,simplices), return
    new vertex and simplex arrays for the Loop subdivided mesh.
    """    
    
    #all edges in the mesh
    edges = set()
    
    for s in simplices:        
        edges.update([frozenset(x) for x in combinations(ravel(s),2)])
    
    edge_index_map = {}    
    for n,e in enumerate(edges):
        edge_index_map[e] = len(vertices) + n
    
    edge_vertices = []
    for e in edges:       
        e0,e1 = sorted(e)
        edge_vertices.append(0.5*(vertices[e0] + vertices[e1]))
    
    new_vertices = concatenate((vertices,array(edge_vertices)))
    new_simplices = [] 
    
    for n,s in enumerate(simplices):
        v0,v1,v2 = ravel(s)
        e01 = edge_index_map[frozenset((v0,v1))]
        e12 = edge_index_map[frozenset((v1,v2))]
        e20 = edge_index_map[frozenset((v2,v0))]
        
        new_simplices.append([v0,e01,e20])
        new_simplices.append([v1,e12,e01])
        new_simplices.append([v2,e20,e12])
        new_simplices.append([e01,e12,e20])
        
    new_simplices = array(new_simplices)
    
    return new_vertices,new_simplices






def triangulate_ncube(vertices,indices):
    n_dims  = ndim(indices) - 1
    n_cubes = indices.shape[0]
    n_verts = vertices.shape[0]
    
    if n_dims <= 1:
        #cube mesh only contains edges
        return vertices,indices





    if n_dims > 2:
        raise NotImplementedError('nCube meshes with n > 2 not supported')

    cell_centers = vertices[indices.reshape(n_cubes,-1)].mean(axis=1)

    n_faces = 2*n_dims*n_cubes

    faces = zeros((n_faces,) + (2,)*(n_dims-1),dtype=indices.dtype)

    for i in range(n_dims):
        s0 = [slice(None,None,None)]*(i+1) + [0] + [slice(None,None,None)]*(n_dims-i-1)
        s1 = [slice(None,None,None)]*(i+1) + [1] + [slice(None,None,None)]*(n_dims-i-1)

        faces[(2*i+0)*n_cubes:(2*i+1)*n_cubes] = indices[s0]
        faces[(2*i+1)*n_cubes:(2*i+2)*n_cubes] = indices[s1]

        #this seems to be the correct pattern
        if (n_dims-1-i) % 2 == n_dims % 2:
            #flip 1
            temp = faces[(2*i+1)*n_cubes:(2*i+2)*n_cubes,0].copy()
            faces[(2*i+1)*n_cubes:(2*i+2)*n_cubes,0] = faces[(2*i+1)*n_cubes:(2*i+2)*n_cubes,1]
            faces[(2*i+1)*n_cubes:(2*i+2)*n_cubes,1] = temp
        else:
            #flip 0
            temp = faces[(2*i+0)*n_cubes:(2*i+1)*n_cubes,0].copy()
            faces[(2*i+0)*n_cubes:(2*i+1)*n_cubes,0] = faces[(2*i+0)*n_cubes:(2*i+1)*n_cubes,1]
            faces[(2*i+0)*n_cubes:(2*i+1)*n_cubes,1] = temp


    face_vertices,face_indices = triangulate_ncube(vertices,faces)

    center_indices = (arange(n_cubes) + face_vertices.shape[0]).reshape((n_cubes,1))
    center_indices = tile(center_indices,(face_indices.shape[0]/n_cubes,1))

    new_vertices = vstack((face_vertices,cell_centers))
    new_indices  = hstack((center_indices,face_indices))

    return new_vertices,new_indices


