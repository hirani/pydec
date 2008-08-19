__all__ = ['is_wellcentered','circumcenter','circumcenter_barycoords']

from numpy import bmat, hstack, vstack, dot, sqrt, ones, zeros, sum, \
        asarray
from numpy.linalg import solve,norm

def is_wellcentered(pts, tol=1e-8):
    """
    Determines whether the M points in N dimensions define a 
    well-centered simplex.
    """
    bary_coords = circumcenter_barycoords(pts)    
    return min(bary_coords) > tol

def circumcenter_barycoords(pts):
    """
    Computes the barycentric coordinates of the circumcenter M, N-dimensional
    points (1 <= M <= N + 1 and N >= 1). The points are given by the rows of 
    an (M)x(N) dimensional matrix pts.
    
    Uses an extension of the method described here:
    http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
    """    

    pts = asarray(pts)

    rows,cols = pts.shape

    assert(rows <= cols + 1)    

    A = bmat( [[ 2*dot(pts,pts.T), ones((rows,1)) ],
               [  ones((1,rows)) ,  zeros((1,1))  ]] )

    b = hstack((sum(pts * pts, axis=1),ones((1))))
    x = solve(A,b)
    bary_coords = x[:-1]  

    return bary_coords
    
def circumcenter(pts):
    """
    Computes the circumcenter and circumradius of M, N-dimensional
    points (1 <= M <= N + 1 and N >= 1). The points are given by the rows of 
    an (M)x(N) dimensional maatrix pts.  
    
    Returns a tuple (center, radius) where center is a
    column vector of length N and radius is a scalar.
        
        In the case of four points in 3D, pts is a 4x3 matrix arranged as:
            
        pts = [ x0 y0 z0 ]
              [ x1 y1 z1 ]
              [ x2 y2 z2 ]          
              [ x3 y3 z3 ]
        
        with return value ([ cx cy cz ], R)    
    
    Uses an extension of the method described here:
    http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
    """
    pts = asarray(pts)      
    bary_coords = circumcenter_barycoords(pts)
    center = dot(bary_coords,pts)
    radius = norm(pts[0,:] - center)
    return (center,radius)
    
    
    
