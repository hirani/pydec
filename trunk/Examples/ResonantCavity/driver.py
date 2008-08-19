"""
Solve the resonant cavity problem with Whitney forms 
and the covolume method.

References:
    Arnold, Douglas N., 
        "Differential Complexes and Numerical Stability"
        Proceedings of the International Congress of Mathematicians, 
        Beijing 2002, Volume 1 : Plenary Lectures, 2002.

"""


from pydec.dec import SimplicialComplex, d, delta
from pydec.fem import whitney_innerproduct
from pydec.io  import read_mesh
from pydec.mesh import loop_subdivision

from scipy import transpose, real, average, zeros
import scipy


#Read in mesh data from file
mesh = read_mesh('mesh_example.xml')

vertices  = mesh.vertices
triangles = mesh.elements


#perform loop subdivision
#vertices,triangles = loop_subdivision(matrix(vertices),matrix(triangles))


sc = SimplicialComplex((vertices,triangles))


#Construct stiffness and mass matrices 
K_fem = transpose(sc[1].d) * whitney_innerproduct(sc,2) * sc[1].d
M_fem = whitney_innerproduct(sc,1)

#Same problem with covolume hodge star
oneform_basis = sc.get_cochain_basis(1)
K_dec = delta(d(oneform_basis)).v


#Eliminate Boundaries from matrices
boundary_edges = sc.boundary()
non_boundary_edges   = set(sc[1].simplex_to_index.keys()) - set(boundary_edges)
non_boundary_indices = [sc[1].simplex_to_index[e] for e in non_boundary_edges]

# eliminate boundary conditions
K_fem = K_fem[non_boundary_indices,:][:,non_boundary_indices]
M_fem = M_fem[non_boundary_indices,:][:,non_boundary_indices]
K_dec = K_dec[non_boundary_indices,:][:,non_boundary_indices]

#TODO consider using sparse eigensolver instead
K_fem = K_fem.todense()
M_fem = M_fem.todense()
K_dec = K_dec.todense()

#Compute eigenvalues and eigenvectors (no sparse eigenvalue solver yet)
fem_eigenvalues,fem_eigenvectors = scipy.linalg.eig(K_fem,M_fem) 
dec_eigenvalues,dec_eigenvectors = scipy.linalg.eig(K_dec)




#Plot the eigenvalues and eigenmodes of the resonant cavity
from pydec.vis import simplex_quivers, triplot
from matplotlib.pylab import plot, quiver, figure, title, show, hold, legend

#number of eigenvalues to plot
NUM_EIGS = 50

fem_values = sorted([x for x in real(fem_eigenvalues) if x > 1e-10])[0:NUM_EIGS]
dec_values = sorted([x for x in real(dec_eigenvalues) if x > 1e-10])[0:NUM_EIGS]


figure()
title('First ' + str(len(fem_values)) + ' Eigenvalues')
hold(True)
plot(fem_values,'ko',label='fem values')
plot(dec_values,'bo',label='dec values')
legend()



#Which non-zero eigenvector to plot?
N = 2


#plot the FEM eigenfield
fem_non_zero_values = real(fem_eigenvectors[:,list(fem_eigenvalues).index(fem_values[N])])
fem_all_values = zeros((sc[1].num_simplices,))
fem_all_values[non_boundary_indices] = fem_non_zero_values
fem_bases,fem_arrows = simplex_quivers(sc,fem_all_values)

figure()
title('FEM Mode #' + str(N+1))
quiver(fem_bases[:,0],fem_bases[:,1],fem_arrows[:,0],fem_arrows[:,1])
triplot(sc.vertices,sc.simplices)



#plot the DEC eigenfield
dec_non_zero_values = real(dec_eigenvectors[:,list(dec_eigenvalues).index(dec_values[N])])
dec_all_values = zeros((sc[1].num_simplices,))
dec_all_values[non_boundary_indices] = dec_non_zero_values
dec_bases,dec_arrows = simplex_quivers(sc,dec_all_values)

figure()
title('DEC Mode #' + str(N+1))
quiver(dec_bases[:,0],dec_bases[:,1],dec_arrows[:,0],dec_arrows[:,1])
triplot(sc.vertices,sc.simplices)

# show the plots
show()

