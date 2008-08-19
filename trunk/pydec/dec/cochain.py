__all__ = ['cochain','Cochain','d','star','delta','laplace_beltrami','laplace_derham']

from scipy import sparse
from pydec.mesh import Simplex

class cochain:
    """
    Represents a cochain associated with a simplical complex
    
    The v member of the cochain is left uninitialized.  This allows functions like
    d(.) and star(.) to operate on single cochains, or groups of cochains together.
    The values associated with each cochain are stored in the columns of v. This is
    especially useful when v is the identity, and thus represents a basis for all cochains.
    Applying operations to this cochain basis allows one to automatically compose the 
    associated matrix operators.
  
    Use the get_cochain() and get_cochain_basis() members of the SimplicialComplex class to
    safely avoid issues with v.    
    """
    
    def __init__(self,complex,dimension,is_primal):                
        self.complex = complex
        self.k = dimension
        self.n = complex.complex_dimension()
        self.is_primal = is_primal
        self.v = None
    def __add__(self,other):
        assert(self.k == other.k and self.complex == other.complex)
        f = cochain(self.complex,self.k,self.is_primal)
        f.v = self.v + other.v
        return f
    def __sub__(self,other):
        assert(self.k == other.k and self.complex == other.complex)        
        f = cochain(self.complex,self.k,self.is_primal)
        f.v = self.v - other.v
        return f
    def __getitem__(self,key):      
        if isinstance(key,Simplex):
            data = self.complex[len(key) - 1]
            index = data.simplex_to_index[key]
            value = self.v.__getitem__(index)
            if key.parity == data.simplex_parity[index]:
                return value
            else:
                return -value            
        else:
            return self.v.__getitem__(key)
    def __setitem__(self,key,value):
        if isinstance(key,Simplex):
            data = self.complex[len(key) - 1]
            index = data.simplex_to_index[key]            
            if key.parity == data.simplex_parity[index]:
                self.v.__setitem__(index,value)
            else:
                self.v.__setitem__(index,-value)            
        else:
            self.v.__setitem__(key,value)
            
    def __str__(self):
        return 'cochain(k='+str(self.k) + ',n=' + str(self.n) + ',is_primal=' + str(self.is_primal) + '\n' + str(self.v) + ')'
    
            
def d(f):
    """
    Implements the discrete exterior derivative d(.)
    
    Accepts a cochain and returns the discrete d applied to the cochain    
    """
    if f.is_primal:
        df = cochain(f.complex, f.k + 1, f.is_primal)
        if f.k == -1:
            df.v = sparse.csr_matrix((f.complex[0].num_simplices,1)) * f.v
        elif f.k < -1 or f.k > f.n + 1:
            df.v = sparse.csr_matrix((1,1)) * f.v
        else:
            df.v = f.complex[f.k].d * f.v
        return df
    else:
        df = cochain(f.complex,f.k + 1,f.is_primal)
        if f.k == -1:
            df.v = sparse.csr_matrix((f.complex[f.n].num_simplices,1)) * f.v
        elif f.k < -1 or f.k > f.n + 1:
            df.v = sparse.csr_matrix((1,1)) * f.v
        else:
            df.v = f.complex[f.n - f.k].boundary * f.v
            df.v *= (-1) ** f.k        
        return df

def star(f):
    """
    Implements the discrete Hodge star *(.)
    
    Accepts a cochain and returns the Hodge star applied to the cochain    
    """
    if f.k == -1 or f.k == f.n + 1:
        starf = cochain(f.complex, f.n - f.k, not f.is_primal)
        starf.v = f.v
        return starf        
    elif f.is_primal:
        starf = cochain(f.complex, f.n - f.k, not f.is_primal)        
        starf.v = f.complex[f.k].star * f.v
        return starf
    else:
        starf = cochain(f.complex, f.n - f.k, not f.is_primal)
        starf.v = f.complex[f.n - f.k].star_inv * f.v
        return starf

def delta(f):
    """
    Implements the discrete codifferental  \delta(.)
    
    Accepts a cochain and returns the codifferental of the cochain    
    """
    sdsf = star(d(star(f)))    
    sdsf.v *= (-1)**(f.n*(f.k-1)+1)
    return sdsf
    
def laplace_derham(f):    
    """
    Implements the discrete Laplace-de Rham \del(.)
    
    Accepts a cochain and returns the Laplace-de Rham of the cochain    
    
    """
    return d(delta(f)) + delta(d(f))
    
def laplace_beltrami(f):
    """
    Implements the discrete Laplace-Beltrami \del(.) = \delta d
    
    Accepts a cochain and returns the Laplace-Beltrami of the cochain    
    
    In the case of 0-forms, the second term of the Laplace-de Rham d(\delta(.)) is 0
    so the Laplace-Beltrami and Laplace-de Rham will be the same.
    """    
    return delta(d(f))
    
    
    
#for backwards compatibility
Cochain = cochain
