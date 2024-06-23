"""
Reads ascii vertex and element files, writes a pydec mesh and displays it
"""

import numpy
from pydec import SimplicialMesh, write_mesh, read_mesh
from matplotlib.pylab import triplot, show

vertices = numpy.loadtxt("v.txt")
elements = numpy.loadtxt("s.txt",dtype='int32') - 1
mymesh = SimplicialMesh(vertices=vertices,indices=elements)
write_mesh("square_8.xml",mymesh,format='basic')
rmesh = read_mesh("square_8.xml")
triplot(rmesh.vertices[:,0], rmesh.vertices[:,1], rmesh.indices)

show()

