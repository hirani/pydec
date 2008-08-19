#from pydec.testing import *
#
#import numpy
#from scipy.sparse import coo_matrix, csr_matrix
#
#from pydec.math.graph import maximal_independent_set
#
#
#class test_maximal_independent_set(TestCase):
#    def is_MIS(self,graph,mis):
#        mis = set(mis)
#        unmarked = set(range(graph.shape[0]))
#        graph = graph.tocoo()
#
#        for i,j in zip(graph.row,graph.col):
#            if i == j:
#                continue #ignore self loops
#            
#            if i in mis and j in mis:
#                return False #not independent
#            
#            if i in mis and j in unmarked:
#                unmarked.remove(j)
#            if j in mis and i in unmarked:
#                unmarked.remove(i)
#
#        return (unmarked == mis) #check maximality
#                
#                    
#        
#    def check_simple(self):
#        """
#        2x2 regular mesh
#        """        
#        A = numpy.matrix([[ 4., -1., -1.,  0.],
#                          [-1.,  4.,  0., -1.],
#                          [-1.,  0.,  4., -1.],
#                          [ 0., -1., -1.,  4.]])
#
#        graph = csr_matrix(A)
#
#        assert_equal(True,self.is_MIS(graph,maximal_independent_set(graph)))
#
#
#    def check_random(self):
#        numpy.random.seed(0)
#
#        def rand_sparse(m,n,nnz_per_row):
#            """
#            Return a sparse csr with a given number of random nonzero entries per row.
#        
#            The actual number of nonzeros may be less than expected due to overwriting.
#            """   
#            nnz_per_row = min(n,nnz_per_row)
#        
#            rows = numpy.arange(m).repeat(nnz_per_row)
#            cols = numpy.random.random_integers(low=0,high=n-1,size=nnz_per_row*m)
#            vals = numpy.random.random_sample(m*nnz_per_row)
#            return coo_matrix((vals,(rows,cols)),(m,n)).tocsr()
#        
#        for n in [2,10,20,50,200]:
#            G = rand_sparse(n,n,min(n/2,10))
#            G = G + G.T
#            assert_equal(True,self.is_MIS(G,maximal_independent_set(G)))
#
