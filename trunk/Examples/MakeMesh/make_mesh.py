"""
Reads ascii vertex and element files, writes a pydec mesh and displays it
"""

import scipy
from pydec import *
from matplotlib.pylab import show

vertices = scipy.loadtxt("v.txt")
elements = scipy.loadtxt("s.txt",dtype='int32') - 1
mymesh = SimplicialMesh(vertices=vertices,indices=elements)
write_mesh("square_8.xml",mymesh,format='basic')
rmesh = read_mesh("square_8.xml")
triplot(rmesh.vertices,rmesh.indices)

show()

