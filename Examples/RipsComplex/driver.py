"""
Detection of holes in coverage in an idealized abstraction of a sensor
network. Hodge decomposition of a random cochain is used to compute a
harmonic cochain.

"""
from scipy import rand
from scipy.linalg import norm
from scipy.sparse.linalg import cg

from pydec import kd_tree, flatten, triplot, read_array
from pydec.dec.rips_complex import rips_complex


pts = read_array('300pts.mtx')  # 300 random points in 2D
rc = rips_complex( pts, 0.15 )

simplices = rc.simplices

cmplx = rc.chain_complex()  # boundary operators [ b0, b1, b2 ]
b1 = cmplx[1].astype(float)  # edge boundary operator
b2 = cmplx[2].astype(float)  # face boundary operator

x = rand(b1.shape[1]) # random 1-chain
# Decompose x using discrete Hodge decomposition
alpha = cg( b1 * b1.T, b1   * x, tol=1e-8)[0]
beta  = cg( b2.T * b2, b2.T * x, tol=1e-8)[0]
h = x - (b1.T * alpha) - (b2 * beta) # harmonic component of x
h /= abs(h).max() # normalize h


#print 'norm( h )',norm(h)
#print 'norm(b1   * h)', norm(b1   * h)  #should be small
#print 'norm(b2.T * h)', norm(b2.T * h)  #should be small


# plot the results
from pylab import figure, title, plot, axis, xlim, ylim, show
from pydec import triplot, lineplot

# Nodes of the Rips complex
figure()
plot(pts[:,0],pts[:,1],'ko')
axis('off')
title('Nodes of the Complex')

# Rips complex (2d triangle mesh)
figure()
title('The Rips Complex')
triplot(pts,simplices[2]) 

# Harmonic 1-chain
figure()
title('Harmonic 1-chain of the Rips complex')
lineplot(pts,simplices[1],linewidths=(15*abs(h)).clip(1,15) ) 

# fiddle with the figures
for i in range(1,4):
    figure(i)
    axis('off')
    axis('equal')
    xlim(-0.1,1.1)
    ylim(-0.1,1.1)

show()

