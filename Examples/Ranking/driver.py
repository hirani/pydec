"""
Least squares ranking on graphs.

Reference :

Least Squares Ranking on Graphs, Hodge Laplacians, Time Optimality,
and Iterative Methods
A. N. Hirani, K. Kalyanaraman, S. Watts
See arXiv:1011.1716v1 [cs.NA] on http://arxiv.org/abs/1011.1716

"""
from numpy import loadtxt
from pydec import abstract_simplicial_complex
from scipy.sparse.linalg import lsqr

# Load graph and edge values
data = loadtxt('data.txt').astype(int)
edges = data[:,:2]

# Create abstract simplicial complex from edges
asc = abstract_simplicial_complex([edges])

omega = data[:,-1] # pairwise comparisons
B1 = asc.chain_complex()[1] # boundary matrix
alpha = lsqr(B1.T, omega)[0] # solve least squares problem

# Set the minimum to 0
alpha = alpha - alpha.min()

print(alpha)






