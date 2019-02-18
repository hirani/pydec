__all__ = ['kd_tree']

from math import sqrt
from heapq import heappush,heappop

class kd_tree:
    class node:
        def point_distance(self,point):
            return sqrt(sum([ (a - b)**2 for (a,b) in zip(point,self.point)]))

        def separator_distance(self,point):
            return point[self.axis] - self.point[self.axis]

    def __repr__(self):
        output = ""
        return "kd_tree< %s points in %s-dimensions >"% (self.num_points,self.k)
    
    def __init__(self, points, values=None):
        """kD-Tree spatial data structure
    
        Parameters
        ----------
        points : array-like
            An N-by-K array of N point coordinates in K dimensions

        Optional Parameters
        -------------------
        values : array-like
            A sequence of N elements associated with the points.
            By default, the integers [0,1,...N-1] are used.

        Examples
        --------
        >>> points = [[0,0],[1,0],[0,1],[1,1]]
        >>> values = ['A','B','C','D']
        >>> kd = kd_tree(points, values)
        >>> kd
        kd_tree< 4 points in 2-dimensions >
        >>> kd.nearest([2,0])
        'B'
        >>> kd.nearest_n([2,0],2)
        ['B', 'D']
        >>> kd.in_sphere([0.1,0.2], 1.1)
        ['A', 'C', 'B']

        """

        lengths = [len(p) for p in points]
        min_dim,max_dim = min(lengths),max(lengths)
        if min_dim != max_dim:
            raise ValueError('points must all have the same dimension')

        if values is None:
            values = range(len(points))
        if len(points) != len(values):
            raise ValueError('points and values must have the same lengths')
        
        self.k = min_dim
        self.num_points = len(points)
        
        self.root = self.__build(zip(points,values),depth=0)

    def __build(self, pv_pairs, depth):
        if not pv_pairs:
            return None
        
        axis = depth % self.k #cycle axis
        
        pv_pairs = sorted(pv_pairs, key=lambda x: x[0][axis])
        
        mid = len(pv_pairs) // 2

        node = self.node()
        node.axis  = axis
        node.point = pv_pairs[mid][0]
        node.value = pv_pairs[mid][1]
        node.left_child  = self.__build(pv_pairs[:mid],   depth+1)
        node.right_child = self.__build(pv_pairs[mid+1:], depth+1)
        return node

    def nearest(self, point, max_dist=float('inf')):
        """Returns the value associated with the nearest points to a given location
        
        Parameters
        ----------
        point : array-like
            Location in space, e.g. [1.5, 2.0]

        Optional Parameters
        -------------------
        max_dist : float
            Ignore points farther than max_dist away from the query point.

        Returns
        -------
        value : single element
            The value associated with the point nearest to the query point.
            Returns None if no points lie within max_dist of the query point
            or the tree is empty.

        """

        x = self.nearest_n(point,n=1,max_dist=max_dist) #list with 0 or 1 elements

        if len(x) == 0:
            return None
        else:
            return x[0]

    def in_sphere(self, point, radius, max_points=None):
        """Returns the values of all points in a given sphere

        Parameters
        ----------
        point : array-like
            Center of the sphere, e.g. [1.5, 2.0]
        radius : float
            Radius of the sphere, e.g. 0.3

        Optional Parameters
        -------------------
        max_points : integer
            An upper-bound on the number of points to return.

        Returns
        -------
        values : list
            List of values associated with all points in the sphere
            defined by point and radius.

        """
        if max_points is None:
            max_points = float('inf')

        return self.nearest_n(point, n=max_points, max_dist=radius)
    

    def nearest_n(self, point, n, max_dist=float('inf')):
        """Returns the values of the nearest n points to a given location
        
        Parameters
        ----------
        point : array-like
            Location in space, e.g. [1.5, 2.0]
        n : integer
            (Maximum) Number of values to return.  Will return
            fewer than n values if the kd_tree contains fewer 
            than n points.

        Optional Parameters
        -------------------
        max_dist : float
            Ignore points farther than max_dist away from the query point.

        Returns
        -------
        values : list
            List of values associated with the n nearest points to
            the query location.

        """

        heap = []
        self.__nearest_n(point, n, max_dist, self.root, heap)
        heap.sort()
        return [ node.value for (neg_dist,node) in reversed(heap) ]

    def __nearest_n(self,point,n,max_dist,current,heap):
        if current is None:
            return max_dist

        pt_dist  = current.point_distance(point)      #distance to this node's point
        sep_dist = current.separator_distance(point)  #signed distance to this node's separating plane

        if pt_dist < max_dist:
            heappush(heap,(-pt_dist,current)) #add this point to the queue
            if len(heap) > n:
                heappop(heap)
            if len(heap) == n:
                max_dist = min(-heap[0][0],max_dist)

        
        if sep_dist < 0:
            max_dist = self.__nearest_n(point,n,max_dist,current.left_child,heap)
        else:
            max_dist = self.__nearest_n(point,n,max_dist,current.right_child,heap)

        if abs(sep_dist) < max_dist:
            #explore other subtree
            if sep_dist < 0:
                return self.__nearest_n(point,n,max_dist,current.right_child,heap)
            else:
                return self.__nearest_n(point,n,max_dist,current.left_child,heap)
        else:
            return max_dist
            
##def inorder(x):
##    if x is not None:
##        return inorder(x.left_child) + [x.value] + inorder(x.right_child)
##    else:
##        return []

