#__all__ = ['matlab_to_complex','complex_to_matlab','sparse_to_ijv']
#
#from pydec.dec import SimplicialComplex,d,star,delta
#from pydec.dec.cochain import Cochain
#
##from scipy import concatenate,random,rand,sparse,zeros,shape,Int,array,matrix,arange,ArrayType
#from scipy import *
#from scipy.io.mio import loadmat,savemat
#
#
#def matlab_to_complex(filename,vertex_array_name = 'v',simplex_array_name='s'):
#    """
#    Load a complex from a MAT file
#    SciPy only supports MAT v5, so save from Matlab with the -V4 option
#    
#        save filename var1 var2 -V4
#    """    
#    
#    dict = {}
#    loadmat(filename,dict)
#    v = dict[vertex_array_name]
#    s = dict[simplex_array_name]
#    
#    s = s.astype(int32)
#    s -= 1
#        
#    for name,arr in dict.iteritems():
#        print name,shape(arr)
#        
#        
#    return SimplicialComplex(v,s)
#
#def complex_to_matlab(filename,complex):
#    """
#    Write a complex and all associated operators to a MAT file
#    """
#    
#    mat_dict = {}
#    
#    ## Export Operators in IJV format
#    for dim in range(complex.complex_dimension + 1):
#        primal_f = complex.get_cochain_basis(dim,True)        
#        dual_f   = complex.get_cochain_basis(dim,False)        
#        
#        mat_dict['primal_d'+str(dim)] = sparse_to_ijv(d(primal_f).v)       
#        mat_dict['dual_d'+str(complex.complex_dimension - dim)] = sparse_to_ijv(d(dual_f).v)
#        
#        mat_dict['primal_star'+str(dim)] = sparse_to_ijv(star(primal_f).v)       
#        mat_dict['dual_star'+str(complex.complex_dimension - dim)] = sparse_to_ijv(star(dual_f).v)
#    
#
#    ##Change 0-based indexing to 1-based
#    for ijv in mat_dict.itervalues():        
#        ijv[:,[0,1]] += 1
#        
#        
#    for dim in range(1,complex.complex_dimension):
#        i_to_s = complex[dim].index_to_simplex
#        s_to_i = complex[dim].simplex_to_index
#        
#        i2s = concatenate([ array([list(i_to_s[x])]) for x in sorted(i_to_s.keys())])
#        i2s += 1
#        mat_dict['sigma'+str(dim)] = i2s.astype(float64)
#    
#    
#    ## 0 and N handled as special cases, 
#    ## 0 because not all verticies are the face of some simplex
#    ## N because the topmost simplices may have an orientation that should be preserved
#    mat_dict['sigma0'] = array(matrix(arange(1,len(complex.vertices)+1)).transpose().astype(float64))           
#    mat_dict['sigma'+str(complex.complex_dimension)] = complex.simplices.astype(float64) + 1
#        
#    mat_dict['v'] = complex.vertices
#    mat_dict['s'] = complex.simplices.astype(float64) + 1
#        
#        
#    
#    
#    savemat(filename,mat_dict)
#
#    
#def sparse_to_ijv(sparse_matrix):
#    """
#    Convert a sparse matrix to a ijv representation.
#    For a matrix with N non-zeros, a N by 3 matrix will be returned
#    
#    Row and Column indices start at 0
#
#    If the row and column entries do not span the matrix dimensions, an additional 
#    zero entry is added for the lower right corner of the matrix
#    """
#    csr_matrix = sparse_matrix.tocsr()
#    ijv = zeros((csr_matrix.size,3))
#    
#    max_row = -1
#    max_col = -1
#    for ii in xrange(csr_matrix.size):
#        ir, ic = csr_matrix.rowcol(ii)
#        data = csr_matrix.getdata(ii)
#        ijv[ii] = (ir,ic,data)    
#        max_row = max(max_row,ir)
#        max_col = max(max_col,ic)
#    
#    
#    rows,cols = shape(csr_matrix)
#    if max_row != (rows - 1) or max_col != (cols - 1):
#        ijv = concatenate((ijv,array([[rows-1,cols-1,0]])))
#        
#    return ijv
#    
#    
#    
#    
#import unittest
#
#class Test_sparse_to_ijv(unittest.TestCase):
#    def setUp(self):
#        random.seed(0) #make tests repeatable  
#        
#    def testsparse_to_ijv(self):
#        cases = []
#        cases.append(((1,1),[(0,0)]))
#        cases.append(((1,3),[(0,0),(0,2)]))
#        cases.append(((7,1),[(5,0),(2,0),(4,0),(6,0)]))
#        cases.append(((5,5),[(0,0),(1,3),(0,4),(0,3),(3,2),(2,0),(4,3)]))
#        
#        for dim,l in cases:
#            s = sparse.lil_matrix(dim)
#            for r,c in l:
#                s[r,c] = 1
#            ijv = sparse_to_ijv(s)
#            
#            self.assertEqual(shape(ijv),(len(l),3))
#            
#            for i,j,v in ijv:
#                self.assert_((i,j) in l)
#                self.assertEqual(v,1)
#                
#            
#class TestFile(unittest.TestCase):
#    def setUp(self):
#        pass
#        
#    def testMatlab(self):
#        sc = matlab_to_complex("../resources/matlab/meshes/unitSqr14")
#        complex_to_matlab("/home/nathan/Desktop/unitSqr14_out",sc)
#        
#        
#    
#if __name__ == '__main__':
#    unittest.main()
