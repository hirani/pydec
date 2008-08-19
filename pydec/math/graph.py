#__all__ = ['maximal_independent_set']
#
#import scipy,numpy
#from numpy.random import permutation
#
#def maximal_independent_set(graph):
#    """
#    Compute a maximal independent set (MIS) of an undirected graph
#    
#        graph:  sparse matrix in CSR or CSC format
#          - may include diagonal entries
#          - may be symmetric (CSR or CSC)
#             or upper triangular (CSR)
#             or lower-triangular (CSC)
#           
#    Returns an array of sorted indices corresponding forming the MIS.
#    """
#    if not (scipy.sparse.isspmatrix_csr(graph) or scipy.sparse.isspmatrix_csc(graph)):
#        raise TypeError,'expected graph in csr_matrix or csc_matrix format'
#    
#    if graph.shape[0] != graph.shape[1]:
#        raise ValueError,'Matrix dimensions must agree'
#        
#    num_nodes = graph.shape[0]
#
#    set_mask = numpy.ones(num_nodes,dtype='byte')  #use 1 for unseen nodes, 2 for MIS nodes, 0 for non-MIS nodes
#
#    indptr,indices = graph.indptr,graph.indices
#
#    
#    code = """
#        #line 32 "graph.py"
#        for(int i = 0; i < num_nodes; i++){
#           if(set_mask(i) == 0){ continue; }
#           
#           const int row_start = indptr(i);
#           const int row_end   = indptr(i+1);           
#           for(int jj = row_start; jj < row_end; jj++){
#              set_mask(indices(jj)) = 0;
#           }
#           set_mask(i) = 2;
#        }
#        """
#        
#    err = scipy.weave.inline(code,
#                             ['indices', 'indptr', 'num_nodes', 'set_mask'],
#                             type_converters = scipy.weave.converters.blitz)
#
#    return set_mask.nonzero()[0]
#
#
#
#
#
#
#def lloyd_cluster(graph,centers,maxiter=None):
#    if not (scipy.sparse.isspmatrix_csr(graph) or scipy.sparse.isspmatrix_csc(graph)):
#        raise TypeError,'expected graph in csr_matrix or csc_matrix form'
#    
#    if graph.shape[0] != graph.shape[1]:
#        raise ValueError,'Matrix dimensions must agree'
#
#    num_nodes = graph.shape[0]
#
#    try:
#        max_dist = numpy.iinfo(graph.dtype).max
#    except:
#        max_dist = numpy.finfo(graph.dtype).max
#
#    #interpret centers argument
#    try:
#        num_centers = int(centers)
#        cluster_centers = numpy.random.permutation(num_nodes)[:num_centers]
#    except:
#        try:
#            cluster_centers = numpy.asarray(centers)
#            num_centers = len(cluster_centers)
#        except:
#            raise ValueError,'unable to parse argument \'centers\''
#
#
#    if num_centers < 1:
#        raise ValueError,'at least one center is required'
#
#
#    if maxiter is None:
#        maxiter = numpy.inf
#        
#
#
#    cluster_radii   = numpy.empty(num_centers,dtype=graph.dtype)
#    
#    distances =  numpy.empty(num_nodes,dtype=graph.dtype)  # distance to nearest cluster center (1) / distance from cluster perimeter (2)
#    current_clusters = -numpy.ones(num_nodes,dtype='i')    # index of nearest cluster center
#    last_clusters    = current_clusters.copy()             # index of nearest cluster center (last iteration)
#
#
#    indptr,indices,data = graph.indptr,graph.indices,graph.data
#    
#    if data.min() < 0:
#        raise ValueError,'negative distances are not permitted'
#
#
#
#    while True:
#        distances[:]               = max_dist  # reset distances
#        distances[cluster_centers] = 0
#
#        last_clusters[:] = current_clusters
#        
#        current_clusters[:]               = -1
#        current_clusters[cluster_centers] = numpy.arange(len(cluster_centers))
#
#
#        code_pass1 = """
#        #line 111 "graph.py"
#        bool any_changes;
#        do{
#           any_changes = false;
#           for(int i = 0; i < num_nodes; i++){
#              for(int jj = indptr(i); jj < indptr(i+1); jj++){
#                 int neighbor = indices(jj);
#                 if(distances(neighbor) + data(jj) < distances(i)){
#                    distances(i)        = distances(neighbor) + data(jj);
#                    current_clusters(i) = current_clusters(neighbor);
#                    any_changes = true;
#                 }
#              }
#           }
#        } while (any_changes);
#        """
#
#        err = scipy.weave.inline(code_pass1,
#                                 ['data', 'indices', 'indptr', 'num_nodes',
#                                  'distances', 'current_clusters'],
#                                 type_converters = scipy.weave.converters.blitz)
#
#        #note: if multiple components exist, there can be some nodes without a distance
#
#
#        distances[:]     =  max_dist      # reset distances
#
#        code_pass2 = """
#        #line 33 "sparse.py"
#
#        //find perimeter nodes
#        for(int i = 0; i < num_nodes; i++){
#           for(int jj = indptr(i); jj < indptr(i+1); jj++){
#              if(current_clusters(indices(jj)) != current_clusters(i)){
#                 distances(i) = 0;          
#              }
#           }
#        }
#        
#        //propagate front inward
#        while(true){
#           bool any_changes = false;
#           for(int i = 0; i < num_nodes; i++){
#              for(int jj = indptr(i); jj < indptr(i+1); jj++){
#                 int neighbor = indices(jj);
#                 if(distances(neighbor) + data(jj) < distances(i)){
#                    distances(i) = distances(neighbor) + data(jj);
#                    any_changes  = true;
#                 }
#              }
#           }
#           if(!any_changes)
#              break;
#        } 
#
#
#        //in the case of ties, use the old centers
#        for(int i = 0; i < num_centers; i++){
#           cluster_radii(i) = distances(cluster_centers(i));
#        }
#
#        for(int i = 0; i < num_nodes; i++){
#           if( cluster_radii(current_clusters(i)) < distances(i) ){
#              cluster_radii(current_clusters(i))   = distances(i);
#              cluster_centers(current_clusters(i)) = i;
#           }
#        }
#        """
#
#        err = scipy.weave.inline(code_pass2,
#                                 ['data', 'indices', 'indptr',
#                                  'num_nodes', 'num_centers',
#                                  'distances', 'current_clusters',
#                                  'cluster_radii', 'cluster_centers'],
#                                 type_converters = scipy.weave.converters.blitz)
#        if (last_clusters == current_clusters).all():
#            break
#
#
#    return distances,current_clusters,cluster_centers
#
#
#
#
#if __name__ == '__main__':
#    from scipy.sparse import spdiags
#    from pydec import triangulate_ncube,cube_grid,simplicial_complex,triplot
#
#    #v,s = simplicial_grid_2d(50)
#
#    v,s = triangulate_ncube(*cube_grid((50,50)))
#
#    sc = simplicial_complex((v,s))
#    G = sc[0].d.T.tocsr() * sc[0].d
#    G = spdiags([G.diagonal()],[0],G.shape[0],G.shape[1],format=G.format) - G
#
#    #random permutation
#    P = scipy.sparse.coo_matrix((numpy.ones(G.shape[0]),(numpy.arange(G.shape[0]),numpy.random.permutation(G.shape[0])))).tocsr()
#    G = P.T.tocsr() * G * P
#    v[P.indices,:] = v.copy()
#    
#    if True:
#        centers = 10 #maximal_independent_set(G)
#        distances,current_clusters,cluster_centers = lloyd_cluster(G,centers)
#
#        from pylab import *
#        #plot(numpy.bincount(distances.astype('i')),'.')
#        #show()
#
#        G_coo = G.tocoo()
#        perimeter_mask = zeros(G.shape[0],dtype='bool')
#        perimeter_mask[G_coo.row[(current_clusters[G_coo.row] != current_clusters[G_coo.col])]] = True
#
#
#        figure()
#        #triplot(v,s)
#        area = 50
#        scatter(v[:,0],v[:,1],area,current_clusters)
#        scatter(v[cluster_centers,0],v[cluster_centers,1],area,'0.0')
#        #scatter(v[perimeter_mask,0],v[perimeter_mask,1],area/10,'1.0')        
#        show()
