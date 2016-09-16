from pydec.testing import *

from scipy import arange, prod, reshape, rand, random, allclose, ndim, zeros
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix

from pydec.io.arrayio import write_array, read_array


#TODO replace with tempfile

filename = '/tmp/pydec_arrayio_testfile.dat'

class TestArrayIO():
    def setUp(self):	
        random.seed(0) #make tests repeatable   
    
    def tearDown(self):
        import os
        os.remove(filename)        

    def test_dense(self):
        sizes = [(2,2),(3,3),(5,1),(1,5)]
        sizes += [(2,2,2),(4,3,2),(1,1,5),(1,5,1),(5,1,1)]
        for dims in sizes:
            mats = [arange(prod(dims)).reshape(dims),rand(*dims)]    
            for A in mats:
                formats = ['binary','ascii']
                if ndim(A) <= 2: formats.append('basic') #use basic when possible
                for format in formats:
                    write_array(filename,A,format=format)
                    
                    B = read_array(filename)
                    assert_almost_equal(A,B,decimal=12)


    def test_sparse(self):
        sizes = [(2,2),(3,3),(1,10),(10,1),(10,10)]
        for dims in sizes:
            base_mats = []
            base_mats.append((rand(*dims) < 0.5)*rand(*dims))  #random matrix with 50% nnz
            base_mats.append(zeros(dims))                      #empty matrix
            base_mats.append(arange(prod(dims)).reshape(dims))

            mats = []
            for base_mat in base_mats:
                mats.append(csr_matrix(base_mat))
                mats.append(csc_matrix(base_mat))
                mats.append(coo_matrix(base_mat))
            
            
            for A in mats:
                formats = ['binary','ascii']                        
                for format in formats:
                    write_array(filename,A,format=format)
                    
                    B = read_array(filename)
                    assert_almost_equal(A.todense(),B.todense(),decimal=12)
                    assert_equal(type(A),type(B))

