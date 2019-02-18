__all__ = ['triplot','lineplot','lineplot2','cube_quivers','simplex_quivers']

try:
    import matplotlib.collections
    import matplotlib.pylab
except ImportError:
    import warnings
    warnings.warn("matplotlib not installed, some loss of functionality will result")


from scipy import rand,asarray,zeros,empty,average
from pydec import barycentric_gradients,combinations,Simplex
import numpy


def triplot(vertices, indices, labels=False):
    """
    Plot a 2D triangle mesh
    """
    
    vertices,indices = asarray(vertices),asarray(indices)

    #3d tensor [triangle index][vertex index][x/y value]
    triangles = vertices[numpy.ravel(indices),:].reshape((indices.shape[0],3,2))
    
    col = matplotlib.collections.PolyCollection(triangles)
    col.set_facecolor('grey')
    col.set_alpha(0.5)
    col.set_linewidth(1)

    #sub =  subplot(111)
    sub = matplotlib.pylab.gca()
    sub.add_collection(col,autolim=True)
    matplotlib.pylab.axis('off')
    sub.autoscale_view()

    if labels:
        barycenters = numpy.average(triangles,axis=1)
        for n,bc in enumerate(barycenters):
            matplotlib.pylab.text(bc[0], bc[1], str(n), {'color' : 'k', 'fontsize' : 8,
                                                         'horizontalalignment' : 'center',
                                                         'verticalalignment' : 'center'
                                                         })

    #matplotlib.pylab.show()


def lineplot2(tails,heads,labels=False,linewidths=1):
    #vertices,indices = asarray(vertices),asarray(indices)
    
    #3d tensor [segment index][vertex index][x/y value]
    #lines = vertices[numpy.ravel(indices),:].reshape((indices.shape[0],2,2))

    data = empty((len(tails),2,2))
    data[:,0,:] = tails
    data[:,1,:] = heads
    col = matplotlib.collections.LineCollection(data)
    col.set_color('k')
    col.set_linewidth(linewidths)

    #sub =  subplot(111)
    sub = matplotlib.pylab.gca()
    sub.add_collection(col,autolim=True)
    matplotlib.pylab.axis('off')
    sub.autoscale_view()

    if labels:
        barycenters = numpy.average(lines,axis=1)
        for n,bc in enumerate(barycenters):
            matplotlib.pylab.text(bc[0], bc[1], str(n), {'color' : 'k', 'fontsize' : 8,
                                                         'horizontalalignment' : 'center',
                                                         'verticalalignment' : 'center'
                                                         })

    #matplotlib.pylab.show()
    

def lineplot(vertices,indices,labels=False,linewidths=1):
    """
    Plot 2D line segments
    """
    vertices,indices = asarray(vertices),asarray(indices)
    
    #3d tensor [segment index][vertex index][x/y value]
    lines = vertices[numpy.ravel(indices),:].reshape((indices.shape[0],2,2))
    
    col = matplotlib.collections.LineCollection(lines)
    col.set_color('k')
    col.set_linewidth(linewidths)

    #sub =  subplot(111)
    sub = matplotlib.pylab.gca()
    sub.add_collection(col,autolim=True)
    matplotlib.pylab.axis('off')
    sub.autoscale_view()

    if labels:
        barycenters = numpy.average(lines,axis=1)
        for n,bc in enumerate(barycenters):
            matplotlib.pylab.text(bc[0], bc[1], str(n), {'color' : 'k', 'fontsize' : 8,
                                        'horizontalalignment' : 'center',
                                        'verticalalignment' : 'center'
                                        })

    #matplotlib.pylab.show()



def cube_quivers(cmplx,vals):
    N = cmplx.complex_dimension()
    quiver_dirs = zeros((cmplx[2].cube_array.shape[0],N))
    
    edge_to_face = cmplx[2].boundary.T.tocsr()
    edge_to_face.data = numpy.abs(edge_to_face.data)

    num_edges = cmplx[1].cube_array.shape[0]
    
    for i in range(N):
        i_edges = (cmplx[1].cube_array[:,-1] == i)
        i_vals  = zeros(num_edges)
        i_vals[i_edges] = vals[i_edges]
        
        quiver_dirs[:,i] = 0.5*(edge_to_face*i_vals)

    quiver_bases = cmplx[2].cube_array[:,:N] + 0.5

    return quiver_bases,quiver_dirs



def simplex_quivers(sc,form):
    """
    Sample a Whitney 1-form at simplex barycenters
    """

    quiver_bases = average(sc.vertices[sc[-1].simplices],axis=1)
    quiver_dirs  = zeros((sc[-1].num_simplices,sc.embedding_dimension()))

    s_to_i = sc[1].simplex_to_index

    for n,s in enumerate(sc[-1].simplices):
        verts = sorted(s)
        
        d_lambda = barycentric_gradients(sc.vertices[verts,:])
        edges   = [Simplex(x) for x in combinations(s,2)]
        indices = [s_to_i[x] for x in edges]
        values  = [form[i] for i in indices]

        for e,v in zip(combinations(range(len(verts)),2),values):
            quiver_dirs[n,:] += v*(d_lambda[e[1]] - d_lambda[e[0]])


        
    quiver_dirs /= (sc.complex_dimension() + 1)

    return quiver_bases,quiver_dirs



##from scipy import *
##from pydec import *
##from pylab import quiver,show

##v = array([[0,0],[1,0],[0,1]])
##s = array([[0,1,2]])

##sc = SimplicialComplex(v,s)

##b,d = simplex_quivers(sc,array([1.0,0.0,0.0]))


##quiver(b[:,0],b[:,1],d[:,0],d[:,1])
##show()





