from numpy import asarray, eye, outer, inner, dot
from numpy.linalg import norm
from scipy.sparse.linalg import cg


def hodge_decomposition(omega):
    """
    For a given k-cochain \alpha there is a unique decomposition
    
    \omega = d(\alpha) + \delta(\beta) (+) h
    
    for k-1 cochain \alpha, k+1 cochain \beta, and harmonic k-cochain h.
    
    This function returns (non-unique) representatives \beta, \gamma, and h 
    which satisfy the equation above.
    
    Example:
        #decompose a random 1-cochain
        sc = SimplicialComplex(...)
        omega = sc.get_cochain(1)
        omega.[:] = rand(*omega.shape)
        (alpha,beta,h) = hodge_decomposition(omega)

    """
    sc = omega.complex    
    
    K = omega.k

    alpha = sc.get_cochain(K - 1)
    beta  = sc.get_cochain(K + 1)    

    #solve for alpha
    A = delta(d(sc.get_cochain_basis(K - 1))).v
    b = delta(omega).v
    alpha.v = cg( A, b, tol=1e-8 )[0]

    #solve for beta
    A = d(delta(sc.get_cochain_basis(K + 1))).v
    b = d(omega).v
    beta.v = cg( A, b, tol=1e-8 )[0]
    
    #solve for h
    h = omega - d(alpha) - delta(beta)    
    
    return (alpha,beta,h)
   

def ortho(A):
    """Separates the harmonic forms stored in the rows of A using a heuristic
    """
    A = asarray(A)

    for i in range(A.shape[0]):
        j = abs(A[i]).argmax()
     
        v = A[:,j].copy()
        if A[i,j] > 0:
            v[i] += norm(v)
        else:
            v[i] -= norm(v)
        Q = eye(A.shape[0]) - 2 * outer(v,v) / inner(v,v)

        A = dot(Q,A)

    return A




import numpy
from scipy import *
from pydec import *

numpy.random.seed(0) # make results consistent

#Read in mesh data from file
mesh = read_mesh('mesh_example.xml')

vertices  = mesh.vertices
triangles = mesh.elements

# remove some triangle from the mesh
triangles = triangles[list(set(range(len(triangles))) - set([30,320,21,198])),:]

# subdivide mesh
#vertices,triangles = loop_subdivision(vertices,triangles)
#vertices,triangles = loop_subdivision(vertices,triangles)

sc = SimplicialComplex((vertices,triangles))

H = []  # harmonic forms

# decompose 4 random 1-cochains
for i in range(4):
    omega = sc.get_cochain(1)
    omega.v[:] = rand(*omega.v.shape)

    (beta,gamma,h) = hodge_decomposition(omega)

    h = h.v
    for v in H:
        h -= inner(v,h) * v
    h /= norm(h)

    H.append(h)


H = ortho(vstack(H))





# plot the results
from pylab import figure, title, quiver, show
from pydec import triplot, simplex_quivers


for n,h in enumerate(H):
    figure()
    title('Harmonic k-form #%d' % n)
    triplot(vertices,triangles) 
    #axis('off')
    
    bases,dirs = simplex_quivers(sc,h)
    quiver(bases[:,0],bases[:,1],dirs[:,0],dirs[:,1]) #,sqrt(numpy.sum(dirs*dirs,axis=1)))

show()

