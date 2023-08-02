from pydec.testing import *

from scipy import random, array, sqrt, sum

from pydec.math.kd_tree import kd_tree


def test_simple():
    pts = [[0.0,0.0],
           [1.0,0.0],
           [1.0,1.0],
           [0.0,1.0],
           [0.5,0.5]]

    kdt = kd_tree(pts)
    
    assert_equal(kdt.nearest([2.0,2.0]), 2)
    assert_equal(kdt.nearest([0.1,0.2]), 0)
    assert_equal(kdt.nearest([0.9,0.0]), 1)
    assert_equal(kdt.nearest([0.4,0.6]), 4)

def test_values():
    pts = [[0.0,0.0],
           [1.0,0.0],
           [1.0,1.0],
           [0.0,1.0],
           [0.5,0.5]]
    vals = ['a','b','c','d','e']
    
    kdt = kd_tree(pts,vals)
    
    assert_equal(kdt.nearest([2.0,2.0]), 'c')
    assert_equal(kdt.nearest([0.1,0.2]), 'a')
    assert_equal(kdt.nearest([0.9,0.0]), 'b')
    assert_equal(kdt.nearest([0.4,0.6]), 'e')

def test_random():
    random.seed(0) #make tests repeatable                 
    for dim in [1,2,3,4]:
        for num_pts in [1,2,5,10,50]:

            pts = (20*rand(num_pts,dim) - 10)
            
            kdt = kd_tree(pts.tolist())

            for sample in (20*rand(5,dim) - 10).tolist(): 
                #5 sample points per test
                distances = sqrt(sum((pts - sample)**2,axis=1))
                sorted_pairs = sorted(zip(distances,range(len(pts))))

                assert_equal(kdt.nearest(sample),sorted_pairs[0][1])
        
                for n in [1,2,5]:
                    assert_equal(kdt.nearest_n(sample,n),[x[1] for x in sorted_pairs[:n]])

