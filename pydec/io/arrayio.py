__all__ = ['read_array','write_array','read_header']

#from scipy.io import read_array,write_array
from scipy import shape,ndim
import numpy
import scipy
import sys


class ArrayIOException(Exception):
    def __init__(self,msg=''): self.msg = msg
    def __str__(self): return self.msg

class FileFormatError(ArrayIOException): pass
      

class ArrayHeader(dict):    
    def tostring(self):        
        self['version'] = '1.0'
        output = str(len(self)) + '\n'
        for key,value in self.items():
            output += key
            output += '='
            output += str(value)
            output += '\n'
        return output


    
#---user functions---#000000#FFFFFF-------------------------------------------------------
def read_array(fid):
    """
    Read an ndarray or sparse matrix from file.
    
    ASCII formats
        basic
        ndarray
        sparse
    BINARY formats
        ndarray
        sparse
        
        Notes:
            ndarray IO makes use to ndarray.tofile() and fromfile()
            
            The following sparse matrix formats are supported:
                csr_matrix
                csc_matrix
                coo_matrix
    """

    if not hasattr(fid, "read"):  fid = open(fid, "rb")
    
    header = read_header(fid)
    
    try: format = header['format']
    except: raise FileFormatError('File format unspecified in file')
    if format not in ['', 'basic', 'ascii', 'binary']: 
        raise FileFormatError('Unknown format: [' + format + ']')    
    
    if format == 'basic':
        return read_basic(fid,header)
    else:
        try:
            array_type = header['type']
        except KeyError:
            raise FileFormatError('Array type unspecified in file: ['+fid.name+']')
        
        if array_type == 'ndarray':
            return read_ndarray(fid,header)    
        elif array_type == 'sparse':
            return read_sparse(fid, header)
        else:
            raise FileFormatError('Unknown array type: [' + array_type + ']')
    


def write_array(fid,A,format='binary'):
    """
    Write an ndarray or sparse matrix to a file
    
        format may be one of ['basic','ascii','binary']
        
        basic
            - Most human readable
            - Only works for arrays of rank 1 and 2
            - Does not work for sparse matrices
        ascii
            - Somewhat human readable
            - Works for ndarrays and sparse matrices
        binary
            - Fastest format
            - Works for ndarrays and sparse matrices
            - Data stored in LittleEndian
    """   
        
    if format not in ['basic', 'ascii', 'binary']: 
        raise ArrayIOException('Unknown format: ['+format+']')
   
    if not hasattr(fid, "read"): fid = open(fid,'wb')
    
    if type(A) is numpy.ndarray:
        A = numpy.ascontiguousarray(A)  #strided arrays break in write
        if format == 'basic':
            if ndim(A) > 2: raise ArrayIOException('basic format only works for rank 1 or 2 arrays')
            write_basic(fid,A)
        else:            
            write_ndarray(fid,A,format)
    elif scipy.sparse.isspmatrix(A):
        if format not in ['ascii', 'binary']: 
            raise ArrayIOException('sparse matrices require ascii or binary format')
        write_sparse(fid,A,format)
    else:
        try:
            A = asarray(A)
            if format == 'basic':
                if ndim(A) > 2: raise ArrayIOException('basic format only works for rank 1 or 2 arrays')
                write_basic(fid,A)
            else:            
                write_ndarray(fid,A,format)
        except:
            raise ArrayIOException('Unknown data type and unable to convert to numpy.ndarray')
        
        
def read_header(fid):
    """
    Read the header of an array file into a dictionary
    """
    if not hasattr(fid, "read"): fid = open(fid, 'rb')
    
    first_line = fid.readline().decode()
    try:    numlines = int(first_line)
    except: 
        print('firstline error: ' + first_line)
        raise ArrayIOException()
    
    #numlines = int(fid.readline())
    header = ArrayHeader()
    for i in range(numlines):
        line = fid.readline().decode().rstrip()
        parts = line.split('=')
        if len(parts) != 2: 
            raise FileFormatError('File header error: line #' + str(i) + ' [' + line + ']')
        header[parts[0]] = parts[1]
    return header
    


#---basic---#000000#FFFFFF------------------------------------------------------
def basic_header(A):
    header = ArrayHeader()   
    header['dims'] = ','.join(list(map(str, A.shape)))
    header['dtype'] = A.dtype.name
    return header

def read_basic(fid,header):
    try:    dimensions = split_on_comma(header['dims'])
    except: raise FileFormatError('Unable to determine dims')

    try: dtype = numpy.typeDict[header['dtype']]
    except: raise FileFormatError('Unable to determine dtype')
    
    if len(dimensions) != 2: raise FileFormatError('basic format only supports 2d arrays')
    if min(dimensions) < 1: raise FileFormatError('all dimensions must be positive')      
   
    return numpy.fromfile(fid,dtype=dtype,count=numpy.prod(dimensions),sep=' ').reshape(dimensions)  
    
def write_basic(fid,A):    
    A = numpy.atleast_2d(A) #force 1d arrays to 2d
    header = basic_header(A)
    header['format'] = 'basic'
    fid.write(header.tostring().encode('utf-8'))
    for row in A:
        row.tofile(fid,sep=' ',format='%.16g')
        fid.write(b'\n')
    

    
#---ndarray---#000000#FFFFFF-------------------------------------------------   
def ndarray_header(A):
    header = ArrayHeader()
    header['type'] = 'ndarray'
    header['rank'] = ndim(A)
    header['dims'] = ','.join(map(str,A.shape))
    header['dtype'] = A.dtype.name
    return header

def read_ndarray(fid,header):
    try:    rank = int(header['rank'])
    except: raise FileFormatError('Unable to determine rank')
    
    try:    dims = split_on_comma(header['dims'])
    except: raise FileFormatError('Unable to determine dims')

    try:    dtype = numpy.typeDict[header['dtype']]
    except: raise FileFormatError('Unable to determine dtype')

    try:    format = header['format']
    except: raise FileFormatError('Unable to determine format')
        
    if len(dims) != rank or min(dims) < 0: 
        raise FileFormatError('Invalid dims')
    
    if format == 'ascii': sep = ' '
    else: sep = ''
   
    if format == 'ascii':
        return numpy.fromfile(fid,dtype=dtype,count=numpy.prod(dims),sep=' ').reshape(dims)
    else:
        A = numpy.fromfile(fid,dtype=dtype,count=numpy.prod(dims),sep='').reshape(dims)
        if sys.byteorder == 'big':
            A = A.byteswap(True) #in-place swap
        return A

def write_ndarray(fid,A,format):
    header = ndarray_header(A)
    header['format'] = format
    fid.write(header.tostring())    

    if format == 'binary':      
        if sys.byteorder == 'little':
            A.tofile(fid)
        else:
            A.byteswap().tofile(fid)
    elif format == 'ascii':        
        A.tofile(fid,sep=' ',format='%.16g')
        if A.size > 0: fid.write('\n') #only introduce newline when something has been written
    else:
        raise ArrayIOException('Unknown file format: ['+format+']')



#---sparse---#000000#FFFFFF-----------------------------------------------------
supported_sparse_formats = ['csr','csc','coo']

def sparse_header(A):
    header = ArrayHeader()
    header['type'] = 'sparse'
    header['sptype'] = A.format
    header['dims'] = ','.join(map(str,A.shape))    
    return header
    
def read_sparse(fid,header):    
    try:    dims = split_on_comma(header['dims'])
    except: raise FileFormatError('Unable to determine dims')

    try:    format = header['sptype']
    except: raise FileFormatError('Unable to determine sparse format')
        
    if len(dims) != 2 or min(dims) < 1: raise FileFormatError('Invalid dims')
    
    if header['sptype'] not in supported_sparse_formats:  
        raise ArrayIOException('Only ' + str(supported_sparse_formats) + ' are supported')
    
    if header['sptype'] == 'csr':
        data   = read_array(fid)
        colind = read_array(fid)
        indptr = read_array(fid)
        return scipy.sparse.csr_matrix((data, colind, indptr), dims)
    elif header['sptype'] == 'csc':
        data   = read_array(fid)
        rowind = read_array(fid)
        indptr = read_array(fid)
        return scipy.sparse.csc_matrix((data, rowind, indptr), dims)
    elif header['sptype'] == 'coo':
        data = read_array(fid)
        row  = read_array(fid)
        col  = read_array(fid)
        return scipy.sparse.coo_matrix((data,(row,col)),dims)

        
def write_sparse(fid,A,format):
    if A.format not in supported_sparse_formats: 
        raise ArrayIOException('Only ' + str(supported_sparse_formats) + ' are supported')
    
    header = sparse_header(A)
    header['format'] = format
    fid.write(header.tostring())
    
    if A.format == 'csr':
        write_array(fid, A.data, format)
        write_array(fid, A.indices, format)
        write_array(fid, A.indptr, format)
    elif A.format == 'csc':
        write_array(fid, A.data,format)
        write_array(fid, A.indices,format)
        write_array(fid, A.indptr,format)
    elif A.format == 'coo':
        write_array(fid, A.data, format)
        write_array(fid, A.row, format)
        write_array(fid, A.col, format)
    else:
        assert(false)

#------Helper functions-------------------------------------------------------------------------
def split_on_comma(to_parse):
    return list(map(int, to_parse.split(',')))
    
