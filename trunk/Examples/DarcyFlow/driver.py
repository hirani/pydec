"""
Darcy flow in a triangle mesh of a planar domain, with constant
velocity prescribed on boundary. 

Reference :

Numerical method for Darcy flow derived using Discrete Exterior Calculus
A. N. Hirani, K. B. Nakshatrala, J. H. Chaudhry
See arXiv:0810.3434v3 [math.NA] on http://arxiv.org/abs/0810.3434

"""
from numpy import mat, zeros, sort, asarray, loadtxt, array, dot, \
     concatenate, sign, vstack, argmax, nonzero
from numpy.linalg import norm, det
from scipy.sparse import bmat
from scipy.sparse.linalg import spsolve
from matplotlib.pylab import figure, gca, triplot
from pydec import simplicial_complex, simplex_quivers, signed_volume

velocity = array([1.,0])
# Read the mesh
vertices = loadtxt('vertices.txt')
triangles = loadtxt('triangles.txt', dtype='int') - 1
# Make a simplicial complex from it
sc = simplicial_complex((vertices,triangles))
# Nk is number of k-simplices
N1 = sc[1].num_simplices
N2 = sc[2].num_simplices
# Permeability is k > 0 and viscosity is mu > 0
k = 1; mu = 1
# The matrix for the full linear system for Darcy in 2D, not taking into
# account the boundary conditions, is :
# [-(mu/k)star1 d1^T ]
# [    d1          Z ] 
# where Z is a zero matrix of size N2 by N2. 
# The block sizes are 
#   N1 X N1    N1 X N2
#   N2 X N1    N2 X N2

d1 = sc[1].d; star1 = sc[1].star
A = bmat([[(-mu/k)*sc[1].star, sc[1].d.T],
          [sc[1].d, None]], format='csr')
b = zeros(N1 + N2) # RHS vector
all_fluxes = zeros(N1)
# Find boundary and internal edges
boundary_edges = sc.boundary(); boundary_edges.sort()
boundary_indices = list(sort([sc[1].simplex_to_index[e] 
                                       for e in boundary_edges]))
num_boundary_edges = len(boundary_indices)
internal_edges = set(sc[1].simplex_to_index.keys()) - set(boundary_edges)
internal_indices = list(sort([sc[1].simplex_to_index[e] 
                                       for e in internal_edges]))
num_internal_edges = sc[1].num_simplices - num_boundary_edges
# Assume triangles oriented the same way so can look at any triangle
s = sign(det(vertices[triangles[0,1:]] - vertices[triangles[0,0]]))
for i, e in enumerate(boundary_edges):
    evector = (-1)**e.parity * (vertices[e[1]] - vertices[e[0]])
    normal = array([-evector[1], evector[0]])
    all_fluxes[boundary_indices[i]] = -s * (1/norm(evector)**2) * \
                              dot(velocity, normal) * \
                              abs(det(vstack((normal, evector))))
pressures = zeros(N2)
# Adjust RHS for known fluxes and pressures
b = b - A * concatenate((all_fluxes,pressures))
# Remove entries of b corresponding to boundary fluxes and known pressure
# Pressure at the right most triangle circumcenter is assumed known
pressure_indices = range(N1, N1+N2)
pressure_indices.remove(N1 + argmax(sc[2].circumcenter[:,0]))
entries_to_keep = concatenate((internal_indices, pressure_indices))
b = b[entries_to_keep]
# Remove corresponding rows and columns of A
A = A[entries_to_keep][:,entries_to_keep]
u = spsolve(A,b)
fluxes = u[0:len(internal_indices)]
pressures[array(pressure_indices)-N1] = u[len(internal_indices):]
# Plot the pressures
figure(); ax = gca();
ax.set_xlabel('x', fontsize=20)
ax.set_ylabel('Pressure', fontsize=20)
ax.plot(sc[2].circumcenter[:,0], pressures, 'ro', markersize=8, mec='r')
ax.hold(True)
# Draw a line of slope -1 through the known presure point
xmax = max(sc[2].circumcenter[:,0])
xmin = min(sc[2].circumcenter[:,0])
ax.plot([xmin, xmax], [xmax - xmin, 0], 'k', linewidth=2)
ax.legend(['DEC', 'Analytical'], numpoints=1)
# Plot the triangles
figure(); ax = gca()
ax.triplot(vertices[:,0], vertices[:,1], triangles)
ax.hold('on')
# Insert the computed fluxes into vector of all fluxes
all_fluxes[internal_indices] = fluxes
# Whitney interpolate flux and sample at barycenters. Then
# rotate 90 degrees in the sense opposite to triangle orientation
v_bases,v_arrows = simplex_quivers(sc, all_fluxes)
v_arrows = -s * vstack((-v_arrows[:,1], v_arrows[:,0])).T
# Plot the resulting vector field at the barycenters
ax.quiver(v_bases[:,0],v_bases[:,1],v_arrows[:,0],v_arrows[:,1],
       units='dots', width=1, scale=1./30)
ax.axis('equal')
ax.set_title('Flux interpolated using Whitney map \n' \
          ' and visualized as velocity at barycenters\n')
