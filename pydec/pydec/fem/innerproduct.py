
__all__ = ['barycentric_gradients','whitney_innerproduct','regular_cube_innerproduct']

from scipy import matrix,sparse,zeros,ones,eye,allclose,alltrue, \
                isreal,real,dot,concatenate,sqrt, \
                arange,array,inner,vstack,atleast_2d,empty,tile, \
                asarray,all,sum,hstack

from scipy.special import factorial, comb
from scipy.linalg import det,inv
from scipy.sparse import coo_matrix
              
from pydec.mesh import simplex
from pydec.math.combinatorial import combinations
from pydec.dec.simplex_array import simplex_array_searchsorted
from pydec.dec.cube_array import cube_array_search,cube_array_boundary

import scipy,numpy


def barycentric_gradients(pts):
    """
    Compute the gradients of the barycentric basis functions over a given simplex
    """            
    V = asarray(pts[1:] - pts[0])
    
    ##all gradients except the first are computed
    grads = dot(inv(inner(V,V)),V) #safer, but slower: grads = scipy.linalg.pinv2(V).T 
    
    ##since sum of all gradients is zero, simply compute the first from the others        
    return vstack((atleast_2d(-numpy.sum(grads,axis=0)),grads))


def massmatrix_rowcols(complex,k):
    """
    Compute the row and column arrays in the COO
    format of the Whitney form mass matrix
    """
    simplices = complex[-1].simplices
    num_simplices = simplices.shape[0]
    p = complex.complex_dimension()
    
    if k == p:
        #top dimension
        rows = arange(num_simplices,dtype=simplices.dtype)
        cols = arange(num_simplices,dtype=simplices.dtype)
        return rows,cols
    
    k_faces = [tuple(x) for x in combinations(range(p+1),k+1)]

    faces_per_simplex = len(k_faces)
    num_faces = num_simplices*faces_per_simplex
    faces     = empty((num_faces,k+1),dtype=simplices.dtype)
   
    for n,face in enumerate(k_faces):
        for m,i in enumerate(face):
            faces[n::faces_per_simplex,m] = simplices[:,i]

    #faces.sort() #we can't assume that the p-simplices are sorted

    indices = simplex_array_searchsorted(complex[k].simplices,faces)

    rows = tile(indices.reshape((-1,1)),(faces_per_simplex,)).flatten()
    cols = tile(indices.reshape((-1,faces_per_simplex)),(faces_per_simplex,)).flatten()

    return rows,cols



def whitney_innerproduct(complex,k):
    """
    For a given SimplicialComplex, compute a matrix representing the 
    innerproduct of Whitney k-forms
    """
    assert(k >= 0 and k <= complex.complex_dimension())    
        
    ## MASS MATRIX COO DATA
    rows,cols = massmatrix_rowcols(complex,k)
    data = empty(rows.shape)


    ## PRECOMPUTATION
    p = complex.complex_dimension()
   
    scale_integration = (factorial(k)**2)/((p + 2)*(p + 1))   
    
    k_forms = [tuple(x) for x in combinations(range(p+1),k)]
    k_faces = [tuple(x) for x in combinations(range(p+1),k+1)]

    num_k_forms = len(k_forms)
    num_k_faces = len(k_faces)

    k_form_pairs = [tuple(x) for x in combinations(k_forms,2)] + [(x,x) for x in k_forms]
    num_k_form_pairs = len(k_form_pairs)
    k_form_pairs_to_index = dict(zip(k_form_pairs,range(num_k_form_pairs)))
    k_form_pairs_to_index.update(zip([x[::-1] for x in k_form_pairs],range(num_k_form_pairs)))
    num_k_face_pairs = num_k_faces**2


    if k > 0:
        k_form_pairs_array = array(k_form_pairs)

    #maps flat vector of determinants to the flattened matrix entries
    dets_to_vals = scipy.sparse.lil_matrix((num_k_face_pairs,num_k_form_pairs))

    k_face_pairs = []
    for face1 in k_faces:
        for face2 in k_faces:
            row_index = len(k_face_pairs)

            k_face_pairs.append((face1,face2))

            for n in range(k+1):
               for m in range(k+1):
                   form1 = face1[:n] + face1[n+1:]
                   form2 = face2[:m] + face2[m+1:]

                   col_index = k_form_pairs_to_index[(form1,form2)]

                   dets_to_vals[row_index,col_index] += (-1)**(n+m)*((face1[n] == face2[m]) + 1)                 

    k_face_pairs_to_index = dict(zip(k_face_pairs,range(num_k_faces**2)))
    dets_to_vals = dets_to_vals.tocsr()
    ## END PRECOMPUTATION


    ## COMPUTATION
    if k == 1:
        Fdet = lambda x : x #det for 1x1 matrices - extend
    else:
        Fdet, = scipy.linalg.flinalg.get_flinalg_funcs(('det',),(complex.vertices,))



    dets = ones(num_k_form_pairs)


    for i,s in enumerate(complex[-1].simplices):

        # for k=1, dets is already correct i.e. dets[:] = 1
        if k > 0:  
            # lambda_i denotes the scalar barycentric basis function of the i-th vertex of this simplex
            # d(lambda_i) is the 1 form (gradient) of the i-th scalar basis function within this simplex
            pts      = complex.vertices[s,:]
            d_lambda = barycentric_gradients(pts)
            
            mtxs = d_lambda[k_form_pairs_array]     # these lines are equivalent to:
            for n,(A,B) in enumerate(mtxs):         #   for n,(form1,form2) in enumerate(k_form_pairs):
                dets[n] = Fdet(inner(A,B))[0]       #       dets[n] = det(dot(d_lambda[form1,:],d_lambda[form2,:].T))

        volume = complex[-1].primal_volume[i]
        vals = dets_to_vals * dets
        vals *= volume * scale_integration  #scale by the volume, barycentric weights, and account for the p! in each whitney form

        #put values into appropriate entries of the COO data array
        data[i*num_k_face_pairs:(i+1)*num_k_face_pairs] = vals


    #now rows,cols,data form a COOrdinate representation for the mass matrix
    shape = (complex[k].num_simplices,complex[k].num_simplices)  
    return coo_matrix((data,(rows,cols)), shape).tocsr()




def regular_cube_innerproduct(rcc,k):      
    """
    For a given regular_cube_complex, compute a matrix
    representing the k-form innerproduct.

    These elements are similar to Whitney forms,
    except using standard linear (bilinear,trilinear,..)
    elements for 0-forms.
    """

    N = rcc.complex_dimension()

    #standard cube is [0,0,..,0] [0,1,...,N]   
    standard_cube  = atleast_2d(array([0]*N + range(N),dtype='i'))
    standard_k_faces = standard_cube
    for i in range(N,k,-1):        
        standard_k_faces = cube_array_boundary(standard_k_faces,i)[0]

        
    k_faces_per_cube = standard_k_faces.shape[0]


    K = zeros((k_faces_per_cube,k_faces_per_cube)) #local stiffness matrix
    h = 1
    V = h**N #cube volume
    scale = V * (1/h)**2 * (1/3.0)**(N-k)
    for i,row_i in enumerate(standard_k_faces):
        for j,row_j in enumerate(standard_k_faces):
            if all(row_i[N:] == row_j[N:]):
                differences = (row_i[:N] != row_j[:N])
                differences[row_i[N:]] = 0                
                K[i,j] = scale * (1.0/2.0)**sum(differences)
            else:
                K[i,j] = 0
        

    CA = rcc[-1].cube_array[:,:N]
    num_cubes = CA.shape[0]

    k_faces  = tile(hstack((CA,zeros((CA.shape[0],k),dtype=CA.dtype))),(1,k_faces_per_cube)).reshape((-1,N+k))
    k_faces += tile(standard_k_faces,(num_cubes,1))
    
    k_face_array = rcc[k].cube_array

    face_indices = cube_array_search(k_face_array,k_faces)

    rows = face_indices.repeat(k_faces_per_cube)
    cols = face_indices.reshape((-1,k_faces_per_cube)).repeat(k_faces_per_cube,axis=0).reshape((-1,))
    data = K.reshape((1,-1)).repeat(num_cubes,axis=0).reshape((-1,))
    
    # temporary memory cost solution - eliminate zeros from COO representation
    nz_mask = data != 0.0
    rows = rows[nz_mask]
    cols = cols[nz_mask]
    data = data[nz_mask]

    shape = (len(k_face_array),len(k_face_array))
    return coo_matrix( (data,(rows,cols)), shape).tocsr()
    
    
  
