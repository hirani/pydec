from pydec.fem.innerproduct import barycentric_gradients

from scipy import zeros,int32,sqrt

        

class ParticleTracer:
    def __init__(self,complex,num_particles):
        
        self.num_particles = num_particles
        self.complex = complex
        self.vertices = complex.vertices
        self.simplices = complex.simplices       
        
        self.particle_positions = zeros((num_particles,complex.embedding_dimension))
        self.particle_velocities = zeros((num_particles,complex.embedding_dimension))       
        
        #the barycenter of each simplex
        self.simplex_barycenters = zeros((len(self.simplices),complex.embedding_dimension))        
        for n,s in enumerate(self.simplices):                
            self.simplex_barycenters[n,:] = self.vertices[sorted(s),:].mean(axis=0)              
        
        #the index of the simplex in which a particle lies
        self.particle_indices = zeros((num_particles),int32)
    
        #fluid velocity within the simplex
        self.simplex_velocites = zeros((len(self.simplices),complex.embedding_dimension))        
    
        #unit normals of the faces of each simplex
        self.face_normals = zeros((len(self.simplices),complex.complex_dimension+1,complex.embedding_dimension))
        
        for n,s in enumerate(self.simplices):
            s_verts = self.vertices[sorted(s),:]            
            bg = -barycentric_gradients(s_verts) #outward face normal is opposite barycentric gradient       
            self.face_normals[n,:,:] = bg/sqrt(sum(multiply(bg,bg),axis=1)) #normalize face normals  
        
        #the first vertex of each face of the simplex
        self.first_face_vertices = zeros((len(self.simplices),complex.complex_dimension+1,complex.embedding_dimension))

        for n,s in enumerate(self.simplices):
            sorted_s = sorted(s)
            self.first_face_vertices[n,:,:] = self.vertices[sorted_s[1:]+[sorted_s[0]],:]               
            
        
        
        #Now compute simplex neighbors, -1 for boundary faces        
        s_to_i = complex[complex.complex_dimension].simplex_to_index
        i_to_s = complex[complex.complex_dimension].index_to_simplex
        f_to_i = complex[complex.complex_dimension - 1].simplex_to_index
        i_to_f = complex[complex.complex_dimension - 1].index_to_simplex
        
        self.simplex_neighbors = zeros((len(s_to_i),complex.complex_dimension+1)) - 1
        
        #initialize to empty sets
        f_to_s = dict([(k,set()) for k in f_to_i.iterkeys()])        
      
        for s in s_to_i.iterkeys():
            for face in simplex_boundary(s):
                f_to_s[face].add(s)        
             
        for s,i in s_to_i.iteritems():
            for n,face in enumerate(simplex_boundary(s)):
                opposite_s = f_to_s[face] - set([s])                
                if len(opposite_s) == 1:
                    self.simplex_neighbors[i,n] = s_to_i[list(opposite_s)[0]]
                assert(len(opposite_s) <= 1)

                    
    def integrate(self, delta_t):        
        assert(delta_t > 0)
        
        #print "PP",self.particle_positions
        
        #fudge factor to ensure that points always lie inside their elements
        EPS = 1e-10
        self.particle_positions = (1.0 - EPS)*self.particle_positions + EPS*self.simplex_barycenters[self.particle_indices]       
        
        for i in range(self.num_particles):
            P = self.particle_positions[i]
            s = self.particle_indices[i]
                        
            time_left = delta_t
            
            while True:                
                V = self.simplex_velocites[s]                        
                N = self.face_normals[s]                    
                FFV = self.first_face_vertices[s]
                N = self.face_normals[s]            
                NV = dot(N,V)
                
            
                #distance to each face                
                D = sum(N*(FFV - P),axis=1)                 
                #speed towards each face
                NV = dot(N,V)    
                NV = select([NV > 0],[NV],0)
                
                CrossingTimes = D/NV
                
                min_face = argmin(CrossingTimes)
                min_time = min(CrossingTimes)  
                
                #print "N",N
                #print "D",D
                #print "NV",NV
                #print "min_time",min_time
                assert(min(D) > 0)
                assert(min_time > 0)
                
                
                if min_time < time_left:
                    #update s,delta_t, position, etc
                    P += min_time * V                    
                    #print "- moving #",str(i)," to ",P                
                    new_s = self.simplex_neighbors[s,min_face]
                    
                    
                    #boundary face, freeze particle 
                    if new_s == -1: 
                        P[:] = (1.0 - EPS)*P + EPS*self.simplex_barycenters[s]
                        break
                    else:
                        s = new_s    
                        P[:] = (1.0 - EPS)*P + EPS*self.simplex_barycenters[s]                        
                        time_left -= min_time
                        
                else:
                    P += time_left * V               
                    self.particle_indices[i] = s
                    assert(allclose(P,self.particle_positions[i]))
                    break
                
from pydec.complex import *
from scipy import *     
import unittest

  
class TestParticleTracer(unittest.TestCase):
    def setUp(self):	
        random.seed(0) #make tests repeatable           

    def test1Din1D(self):
        """
        Flow to the right on the line segment [0,8]
        """
        sc = SimplicialComplex(matrix([[0],[1],[2],[4],[8]]),matrix([[0,1],[1,2],[2,3],[3,4]]))
        pt = ParticleTracer(sc,2)
        
        pt.particle_indices[0] = 0
        pt.particle_indices[1] = 1
        
        pt.particle_positions[0] = 0
        pt.particle_positions[1] = 1.5
        
        pt.simplex_velocites[0] = [1.0]
        pt.simplex_velocites[1] = [1.0]
        pt.simplex_velocites[2] = [1.0]
        pt.simplex_velocites[3] = [1.0]
        
        pt.integrate(0.5)
        assert(allclose(pt.particle_positions,array([[0.5],[2.0]])))                
        pt.integrate(0.25)        
        assert(allclose(pt.particle_positions,array([[0.75],[2.25]])))        
        pt.integrate(1.1)
        assert(allclose(pt.particle_positions,array([[1.85],[3.35]])))        
        pt.integrate(1.0)
        assert(allclose(pt.particle_positions,array([[2.85],[4.35]])))                   
        pt.integrate(4.0)
        assert(allclose(pt.particle_positions,array([[6.85],[8.0]])))
        pt.integrate(2.0)
        assert(allclose(pt.particle_positions,array([[8.0],[8.0]])))
        pt.integrate(1.0)
        assert(allclose(pt.particle_positions,array([[8.0],[8.0]])))
        
        
    def test1Din2D(self):
        """
        Counter-clockwise flow on the border of a unit square
        """
        sc = SimplicialComplex(matrix([[0,0],[1,0],[1,1],[0,1]]),matrix([[0,1],[1,2],[2,3],[3,0]]))
        pt = ParticleTracer(sc,1)
        
        pt.particle_indices[0] = 0      
        
        pt.particle_positions[0] = [0, 0]        
        
        pt.simplex_velocites[0] = [1.0,0.0]
        pt.simplex_velocites[1] = [0.0,1.0]
        pt.simplex_velocites[2] = [-1.0,0.0]
        pt.simplex_velocites[3] = [0.0,-1.0]
        
        pt.integrate(0.5)
        assert(allclose(pt.particle_positions,array([[0.5,0.0]])))                
        pt.integrate(0.6)        
        assert(allclose(pt.particle_positions,array([[1.0,0.1]])))                
        pt.integrate(0.7)
        assert(allclose(pt.particle_positions,array([[1.0,0.8]])))        
        pt.integrate(0.8)
        assert(allclose(pt.particle_positions,array([[0.4,1.0]])))        
        pt.integrate(0.9)
        assert(allclose(pt.particle_positions,array([[0.0,0.5]])))        
        pt.integrate(1.0)
        assert(allclose(pt.particle_positions,array([[0.5,0.0]])))        
        pt.integrate(1.0)
        assert(allclose(pt.particle_positions,array([[1.0,0.5]])))        
        pt.integrate(1.0)
        assert(allclose(pt.particle_positions,array([[0.5,1.0]])))        
        pt.integrate(1.0)
        assert(allclose(pt.particle_positions,array([[0.0,0.5]])))        
        pt.integrate(1.0)
        assert(allclose(pt.particle_positions,array([[0.5,0.0]])))        


    def test2Din2D(self):
        """
        A backwards 'C' shaped flow in a square with side length 2
        
                <-<-<-
                     ^
                     ^
                ->->->        
        """       
        sc = SimplicialComplex(matrix([[0,0],[1,0],[2,0],[0,1],[1,1],[2,1],[0,2],[1,2],[2,2]]), \
                            matrix([[0,1,3],[1,4,3],[1,2,4],[2,5,4],[3,4,6],[4,7,6],[4,8,7],[4,5,8]]))
        pt = ParticleTracer(sc,4)
        
        pt.particle_indices[0] = 0
        pt.particle_indices[1] = 1
        pt.particle_indices[2] = 2
        pt.particle_indices[3] = 3
        
        pt.particle_positions[0] = [1.0/3.0, 1.0/3.0]
        pt.particle_positions[1] = [2.0/3.0, 2.0/3.0]
        pt.particle_positions[2] = [4.0/3.0, 1.0/3.0]
        pt.particle_positions[3] = [5.0/3.0, 2.0/3.0]
                
        pt.simplex_velocites[0] = [1.0, 0.0]
        pt.simplex_velocites[1] = [1.0, 0.0]
        pt.simplex_velocites[2] = [1.0, 0.0]
        pt.simplex_velocites[3] = [0.0, 1.0]
        pt.simplex_velocites[4] = [-1.0, 0.0]
        pt.simplex_velocites[5] = [-1.0, 0.0]
        pt.simplex_velocites[6] = [-1.0, 0.0]
        pt.simplex_velocites[7] = [0.0, 1.0]
        
        
        pt.integrate(1.0/3.0)                
        assert(allclose(pt.particle_positions,array([[2.0/3.0,1.0/3.0],[1.0,2.0/3.0],[5.0/3.0,1.0/3.0],[5.0/3.0,1.0]])))     
        pt.integrate(1.0/3.0)        
        assert(allclose(pt.particle_positions,array([[1.0,1.0/3.0],[4.0/3.0,2.0/3.0],[5.0/3.0,2.0/3.0],[5.0/3.0,4.0/3.0]])))
        pt.integrate(1.0/3.0)        
        assert(allclose(pt.particle_positions,array([[4.0/3.0,1.0/3.0],[4.0/3.0,1.0],[5.0/3.0,1.0],[5.0/3.0,5.0/3.0]])))
        pt.integrate(1.0/3.0)        
        assert(allclose(pt.particle_positions,array([[5.0/3.0,1.0/3.0],[4.0/3.0,4.0/3.0],[5.0/3.0,4.0/3.0],[4.0/3.0,5.0/3.0]])))
        pt.integrate(1.0/3.0)        
        assert(allclose(pt.particle_positions,array([[5.0/3.0,2.0/3.0],[1.0,4.0/3.0],[5.0/3.0,5.0/3.0],[1.0,5.0/3.0]])))
        pt.integrate(1.0/3.0)        
        assert(allclose(pt.particle_positions,array([[5.0/3.0,1.0],[2.0/3.0,4.0/3.0],[4.0/3.0,5.0/3.0],[2.0/3.0,5.0/3.0]])))
        pt.integrate(1.0/3.0)        
        assert(allclose(pt.particle_positions,array([[5.0/3.0,4.0/3.0],[1.0/3.0,4.0/3.0],[1.0,5.0/3.0],[1.0/3.0,5.0/3.0]])))
        pt.integrate(1.0/3.0)        
        assert(allclose(pt.particle_positions,array([[5.0/3.0,5.0/3.0],[0.0/3.0,4.0/3.0],[2.0/3.0,5.0/3.0],[0.0,5.0/3.0]])))
        pt.integrate(1.0/3.0)        
        assert(allclose(pt.particle_positions,array([[4.0/3.0,5.0/3.0],[0.0/3.0,4.0/3.0],[1.0/3.0,5.0/3.0],[0.0,5.0/3.0]])))
        pt.integrate(1.0/3.0)        
        assert(allclose(pt.particle_positions,array([[1.0,5.0/3.0],[0.0/3.0,4.0/3.0],[0.0,5.0/3.0],[0.0,5.0/3.0]])))
        pt.integrate(1.0/3.0)        
        assert(allclose(pt.particle_positions,array([[2.0/3.0,5.0/3.0],[0.0/3.0,4.0/3.0],[0.0,5.0/3.0],[0.0,5.0/3.0]])))
        pt.integrate(1.0/3.0)        
        assert(allclose(pt.particle_positions,array([[1.0/3.0,5.0/3.0],[0.0/3.0,4.0/3.0],[0.0,5.0/3.0],[0.0,5.0/3.0]])))
        pt.integrate(1.0/3.0)        
        assert(allclose(pt.particle_positions,array([[0.0/3.0,5.0/3.0],[0.0/3.0,4.0/3.0],[0.0,5.0/3.0],[0.0,5.0/3.0]])))
        pt.integrate(1.0/3.0)        
        assert(allclose(pt.particle_positions,array([[0.0/3.0,5.0/3.0],[0.0/3.0,4.0/3.0],[0.0,5.0/3.0],[0.0,5.0/3.0]])))


        
        
if __name__ == '__main__':
    unittest.main()
