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
    
    def __init__(self,points,values=None):
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

    def __build(self,pv_pairs,depth):
        if not pv_pairs:
            return None
        
        axis = depth % self.k #cycle axis
        
        pv_pairs.sort(key=lambda x: x[0][axis])
        
        mid = len(pv_pairs) / 2

        node = self.node()
        node.axis  = axis
        node.point = pv_pairs[mid][0]
        node.value = pv_pairs[mid][1]
        node.left_child  = self.__build(pv_pairs[:mid],   depth+1)
        node.right_child = self.__build(pv_pairs[mid+1:], depth+1)
        return node

    def nearest(self,point,max_dist=float('inf')):
        x = self.nearest_n(point,n=1,max_dist=max_dist) #list with 0 or 1 elements

        if len(x) == 0:
            return None
        else:
            return x[0]

    def in_sphere(self,point,radius,max_points=None):
        if max_points is None:
            max_points = float('inf')

        return self.nearest_n(point,n=max_points,max_dist=radius)
    

    def nearest_n(self,point,n,max_dist=float('inf')):
        heap = []
        self.__nearest_n(point,n,max_dist,self.root,heap)
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

