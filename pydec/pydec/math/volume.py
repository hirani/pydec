__all__ = ['unsigned_volume','signed_volume']

from scipy import sqrt,inner,shape,asarray
from scipy.special import factorial
from scipy.linalg import det


def unsigned_volume(pts):
    """Unsigned volume of a simplex    
    
    Computes the unsigned volume of an M-simplex embedded in N-dimensional 
    space. The points are stored row-wise in an array with shape (M+1,N).
    
    Parameters
    ----------
    pts : array
        Array with shape (M+1,N) containing the coordinates
        of the (M+1) vertices of the M-simplex.

    Returns
    -------
    volume : scalar
        Unsigned volume of the simplex

    Notes
    -----
    Zero-dimensional simplices (points) are assigned unit volumes.
        

    Examples
    --------
    >>> # 0-simplex point 
    >>> unsigned_volume( [[0,0]] )
    1.0
    >>> # 1-simplex line segment
    >>> unsigned_volume( [[0,0],[1,0]] )             
    1.0
    >>> # 2-simplex triangle 
    >>> unsigned_volume( [[0,0,0],[0,1,0],[1,0,0]] ) 
    0.5


    References
    ----------
    [1] http://www.math.niu.edu/~rusin/known-math/97/volumes.polyh

    """       
    
    pts = asarray(pts)
    
    M,N = pts.shape
    M -= 1

    if M < 0 or M > N:
        raise ValueError('array has invalid shape')
    
    if M == 0:
        return 1.0 
        
    A = pts[1:] - pts[0]
    return sqrt(abs(det(inner(A,A))))/factorial(M)
    
    
def signed_volume(pts):
    """Signed volume of a simplex    
    
    Computes the signed volume of an M-simplex embedded in M-dimensional 
    space. The points are stored row-wise in an array with shape (M+1,M).
    
    Parameters
    ----------
    pts : array
        Array with shape (M+1,M) containing the coordinates
        of the (M+1) vertices of the M-simplex.

    Returns
    -------
    volume : scalar
        Signed volume of the simplex


    Examples
    --------
    >>> # 1-simplex line segment
    >>> signed_volume( [[0],[1]] )           
    1.0
    >>> # 2-simplex triangle 
    >>> signed_volume( [[0,0],[1,0],[0,1]] ) 
    0.5
    >>> # 3-simplex tetrahedron
    >>> signed_volume( [[0,0,0],[3,0,0],[0,1,0],[0,0,1]] ) 
    0.5

    References
    ----------
    [1] http://www.math.niu.edu/~rusin/known-math/97/volumes.polyh

    """       
    
    pts = asarray(pts)
    
    M,N = pts.shape
    M -= 1

    if M != N:
        raise ValueError('array has invalid shape')
        
    A = pts[1:] - pts[0]
    return det(A)/factorial(M)
