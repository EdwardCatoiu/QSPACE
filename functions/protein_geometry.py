import math
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
from scipy.optimize import leastsq

    

def distance_two_points(p1, p0):
    return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2 + (p0[2] - p1[2])**2)
    
    

def project_q_into_plane(q, n_plane):
    #plane solved from above
    n_plane =  np.array(n_plane)
   
    #find p, point in the plane
    pz_plane = -n_plane[3] / n_plane [2]
    p = np.array([0 , 0 , pz_plane]) #point in the plane (0,0, -D / C)
    #normal vector
    n = np.array([n_plane[0], n_plane[1],n_plane[2]]) #normal vector
    #unit normal vector
    n_norm =  np.array(find_unit_norm(n)) #unit normal vector

    #point of interest
    q =  np.array(q)
    #projection of q into plane
    q_proj = q - np.dot(q-p, n_norm )*n_norm
    #print 'q_proj : ' , q_proj
    #print 'Plane :   %fX + %fY + %fZ + %f = 0' %(n_plane[0],n_plane[1],n_plane[2],n_plane[3])
    return q_proj, n_norm

def find_unit_norm(n):
    #print n
    A = n[0]
    B = n[1]
    C = n[2]
    magnitude = find_magnitude_of_vector(n)
    norm_A =  A / magnitude
    norm_B =  B / magnitude
    norm_C =  C / magnitude
    return np.array([norm_A, norm_B, norm_C])


def find_magnitude_of_vector(n):
    A = n[0]
    B = n[1]
    C = n[2]
    mag = math.sqrt(A**2 + B **2 + C** 2)
    return mag

def find_distance_q_to_n_plane(q, n_plane):
    return np.abs(np.dot(n_plane[0:3], q) +n_plane[3])/find_magnitude_of_vector(n_plane[0:3])

def find_lstsq_plane_equation(x_list,y_list,z_list):
    # coordinates (XYZ) of C1, C2, C4 and C5
    XYZ = np.array([ x_list, y_list,z_list])
    # Inital guess of the plane
    p0 = [10, 1, 10, 1]
    def f_min(X,p):
        plane_xyz = p[0:3]
        distance = (plane_xyz*X.T).sum(axis=1) + p[3]
        return distance / np.linalg.norm(plane_xyz)
    def residuals(params, signal, X):
        return f_min(X, params)
    sol = leastsq(residuals, p0, args=(None, XYZ))[0]
    #print "Solution: ", sol
    #print "Old Error: ", (f_min(XYZ, p0)**2).sum()
    new_error = (f_min(XYZ, sol)**2).sum()
    return sol, new_error

def lstsq_determine_membrane_surface_plane(pdb_coords_per_gene_leaf,
                                           variance_threshold = 0.05 ):
    #input : PDB residue : coords --> for each leaflet!!!
    residue_index_list = list(pdb_coords_per_gene_leaf.keys())


    if len(residue_index_list) < 3:
        return [], []
    #manually find plane of 3 points
    if len(residue_index_list) == 3:
        vectors = list(pdb_coords_per_gene_leaf.values())
        p1 = np.array(vectors[0])
        p2 = np.array(vectors[1])
        p3 = np.array(vectors[2])
        # These two vectors are in the plane
        v1 = p3 - p1
        v2 = p2 - p1
        # the cross product is a vector normal to the plane
        cp = np.cross(v1, v2)
        a, b, c = cp
        # This evaluates a * x3 + b * y3 + c * z3 which equals d
        d = np.dot(cp, p3)    
        n_plane = [a,b,c,d]
        return n_plane, residue_index_list
    #use optimal plane solver    
    if len(residue_index_list) == 4:
        x_list = []
        y_list = []
        z_list = []
        for residue, residue_xyz in pdb_coords_per_gene_leaf.items():
            #residue sweep is the residue we are removing
            x_list.append(residue_xyz[0])
            y_list.append(residue_xyz[1])
            z_list.append(residue_xyz[2])
        #find the equation of the plane and the error associated with the calculation
        sol, new_error = find_lstsq_plane_equation(x_list,y_list,z_list)
        n_plane = sol
        return n_plane, residue_index_list

    #iterative loop over all the residues in the leaf, residue list gets shorter after the
    #most error-prone residue is removed
    area_list = {}
    variance_list = []
    mutual_info_reached = False
    while len(residue_index_list) > 4 :
        res_error = {}
        error_list= []
        for residue_sweep in residue_index_list:
            x_list = []
            y_list = []
            z_list = []
            list_dict = {}

            for residue, residue_xyz in pdb_coords_per_gene_leaf.items():
                #residue sweep is the residue we are removing
                if residue == residue_sweep:
                    continue
                if residue not in residue_index_list:
                    #pass because the residue has been removed from the calculation in a previous loop
                    continue
                x_list.append(residue_xyz[0])
                y_list.append(residue_xyz[1])
                z_list.append(residue_xyz[2])

            #find the equation of the plane and the error associated with the calculation
            sol, new_error = find_lstsq_plane_equation(x_list,y_list,z_list)
            list_dict.update({'error': new_error})
            list_dict.update({'plane_solution': sol})


            #res error holds the information for each residue that has been removed:
            #resulting plane equation and error associated with it 
            res_error.update({residue_sweep : list_dict})
            #error list contains the error resulting from removal of each residue
            #we want to minimize the error (if we remove the residue, the remaining plane has the least error!)
            error_list.append(new_error)

        #we remove the residue from the residue index list if it results in lowest error
        for residue in res_error.keys():
            if res_error[residue].get('error') == min(error_list):
                residue_index_list.remove(residue) #removed

                #here we identify the PLANE EQUATION {AX, BY, CZ, D} of that minimum error residue
                n_plane = res_error[residue]['plane_solution']

        #we now project find all the distances from RESIDUES INTO THE NEW CALCULATED PLANE
        distance_list = []
        for residue, q in pdb_coords_per_gene_leaf.items():
            if residue in res_error.keys():
                distance = find_distance_q_to_n_plane(q, n_plane)
                distance_list.append(distance)

        variance_list.append(np.var(distance_list))
        v = 0
        if len(variance_list) > 1:
            while  v < len(variance_list)-1:
                variance_change =  (variance_list[v] - variance_list[v+1])/variance_list[0] 
                if variance_change < variance_threshold: #manual input threshold
                    mutual_info_reached = True
                    break
                v = v + 1 
        #this means that removal of one more residue will not add enough information to the new plane!
        if  mutual_info_reached:
            break


        #now we must USE PCA to reduce dimentionality, followed by convexHull to outline polygon

    #vert_coords = find_polygon_vertices(q_proj_dict)
    #surface_area = area(vert_coords)
    return n_plane , residue_index_list

def find_angle_between_vector(vector1, vector2):
    return math.degrees(np.arccos(np.abs(np.dot(vector1 , vector2))))
def unit_normal(a, b, c):
    x = det([[1,a[1],a[2]],
             [1,b[1],b[2]],
             [1,c[1],c[2]]])
    y = det([[a[0],1,a[2]],
             [b[0],1,b[2]],
             [c[0],1,c[2]]])
    z = det([[a[0],a[1],1],
             [b[0],b[1],1],
             [c[0],c[1],1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

def det(a):
    return a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][2]*a[1][1]*a[2][0] - a[0][1]*a[1][0]*a[2][2] - a[0][0]*a[1][2]*a[2][1]

def closest_leaf_embeddedness(dfpdb, leaf1,leaf2,vectors):

    membrane_thickness = {}
    closest_dict = {}
    embedded_dict = {}
    for index in (dfpdb.index.tolist()):
        row = dfpdb.loc[index]
        x = row.get('x')
        y = row.get('y')
        z = row.get('z')

        dist_leaf1 = find_distance_q_to_n_plane(q = [x,y,z],n_plane=vectors[leaf1])
        dist_leaf2 = find_distance_q_to_n_plane(q = [x,y,z],n_plane=vectors[leaf2])
        distances = {leaf1 : dist_leaf1, leaf2 : dist_leaf2}
        closest_leaf, furthest_leaf = get_closest_leaf(distances)


        proj_res_leaf1, N_unit_vector  = project_q_into_plane(q = [x,y,z],n_plane=vectors[leaf1])
        distance_between_membranes = find_distance_q_to_n_plane(q = proj_res_leaf1,n_plane=vectors[leaf2])
#         membrane_thickness += [distance_between_membranes]

        if np.max(list(distances.values())) > distance_between_membranes:
            embedded = False
        else:
            embedded = True

        closest_dict.update({index: closest_leaf})
        embedded_dict.update({index: embedded})
        membrane_thickness.update({index: distance_between_membranes})
        
    dfpdb['closest_leaf' ] = pd.Series(closest_dict)
    dfpdb['embedded' ] = pd.Series(embedded_dict)
    dfpdb['membrane_thickness' ] = pd.Series(membrane_thickness)
    return dfpdb

def get_closest_leaf(distance_dict):
    df= pd.DataFrame()
    df['distance']  = pd.Series(distance_dict)
    df = df.sort_values(by = 'distance')
    return df.first_valid_index(), df.last_valid_index()

def get_projected_polygons_to_leaf(dfpdb, vectors, only_use_embedded = True):

   

    
    if only_use_embedded:
        dfpdb = dfpdb[dfpdb.embedded == True]

    
    
    leaf_polygons = {}
    
    for leaf in dfpdb.closest_leaf.unique():
        if leaf not in vectors:
            continue
        dfpdb_leaf = dfpdb[dfpdb.closest_leaf == leaf]
        vector = vectors[leaf]
        projected_residues_dict = {}

        for index, row in dfpdb_leaf.iterrows():
            x = row.get('x')
            y = row.get('y')
            z = row.get('z')

            dist_leaf1 = find_distance_q_to_n_plane(q = [x,y,z],n_plane=vectors[leaf])

            projected_residue, unit_vector  = project_q_into_plane(q = [x,y,z],n_plane=vector)    
            projected_residues_dict.update({index : projected_residue})
        leaf_polygons.update({leaf : projected_residues_dict})
    return leaf_polygons

def get_polygon_vertices(polygon_coords):
    
    model = PCA(n_components= 2).fit(polygon_coords)
    proj_vertices = model.transform(polygon_coords)
    hull_kinda = ConvexHull(proj_vertices)
    vertices =  hull_kinda.vertices
    orig_cords =  model.inverse_transform(proj_vertices)


    polygon_vert_coords = []
    for index_value  in vertices:
        polygon_vert_coords.append(orig_cords[index_value])
    return polygon_vert_coords

def find_area(poly):
    if len(poly) < 3: # not a plane - no area
        return 0

    total = [0, 0, 0]
    for i in range(len(poly)):
        vi1 = poly[i]
        if i is len(poly)-1:
            vi2 = poly[0]
        else:
            vi2 = poly[i+1]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result/2)

def get_closest_leaf_NEW(distance_dict):
    df= pd.DataFrame()
    df['distance']  = pd.Series(distance_dict)
    df = df.sort_values(by = 'distance')
    return df


def closest_leaf_embeddedness_NEW(dfpdb, vectors):

    membrane_thickness = {}
    closest_dict = {}
    embedded_dict = {}
    
    
    
    
    for index in (dfpdb.index.tolist()):
        row = dfpdb.loc[index]
        x = row.get('x')
        y = row.get('y')
        z = row.get('z')

        
        distances = {}
        for leaf, vector in vectors.items():
            
            dist = find_distance_q_to_n_plane(q = [x,y,z],n_plane=vector)
            distances.update({leaf : dist})
        
        dfclosest = get_closest_leaf_NEW(distances)
        closest_leaf = dfclosest.first_valid_index()
        second_leaf = dfclosest.index.tolist()[1]
        

        proj_res_leaf1, N_unit_vector  = project_q_into_plane(q = [x,y,z],n_plane=vectors[closest_leaf])
        distance_between_membranes = find_distance_q_to_n_plane(q = proj_res_leaf1,n_plane=vectors[second_leaf])
#         membrane_thickness += [distance_between_membranes]

        if np.max(list(distances.values())) > distance_between_membranes:
            embedded = False
        else:
            embedded = True

        closest_dict.update({index: closest_leaf})
        embedded_dict.update({index: embedded})
        membrane_thickness.update({index: distance_between_membranes})
        
    dfpdb['closest_leaf' ] = pd.Series(closest_dict)
    dfpdb['embedded' ] = pd.Series(embedded_dict)
    dfpdb['membrane_thickness' ] = pd.Series(membrane_thickness)
    return dfpdb