# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)

from .OPMserver import *

from .membraneMetrics import *



def run_005D1_opmPlanes(potential_membrane_structures,
                      OPMstructureFolder = qspaceDirs['opmOutputStructuresDir'] ,
                       force_rerun = False,
                       OPMJsonFolder = qspaceDirs['opmOutputDataDir'],
                       OPMcsvFolder = qspaceDirs['opmCsvDir'],
                       checked = []
                      ):
    
    planes_dict= {}
    previousData_planes_dict= {}
    json_planes_path = op.join(qspaceDirs['DataOutput_dir'], '005D-MembranePlanesOPM.json')
    if not op.exists(json_planes_path):
        with open(json_planes_path ,'w') as f:
            json.dump(previousData_planes_dict, f)
    
    if op.exists(json_planes_path) and not force_rerun:
        with open(json_planes_path ,'r') as f:
            previousData_planes_dict = json.load(f)
    
    
    print ("Getting OPM membrane data ...\n------------------------------------------------")
    print ("Num AAs \t structureId \t N-leaf AAs \t O-leaf AAs \t Label")

    for s_id in tqdm(potential_membrane_structures):
    
        sfile = op.join(OPMstructureFolder, "OPMstructure-{}out.pdb".format(s_id.replace('-','_')))
        if not op.exists(sfile):
            print ("{} OPM file DNE >> {}".format(s_id, sfile))
            continue   

        with open(json_planes_path ,'r') as f:
            previousData_planes_dict = json.load(f)

        if s_id in previousData_planes_dict and not force_rerun:
            planes_dict.update({s_id: previousData_planes_dict[s_id]})
            print ('Using previous data for .... {}'.format(s_id))
            continue
        dfpdb = get_pdb_coordinates(s_id, sfile)

        closest_dict = {}
        embedded_dict = {}

        N_coords, O_coords = get_coordinates_of_opm_membrane(sfile, num_residues_to_consider  = 50)


        if O_coords == {} or N_coords == {}:  #the OPM pdb does not have membrane atoms attached for 1 (or 2) membranes
    #         print  len(dfpdb) , '\t', s_id, len(N_coords),  len(O_coords)

            if O_coords == {} and N_coords == {}:
                #we need to check the HTML data for transmembrane residues in this case
                outFileJSON = op.join(OPMJsonFolder, "OPMjson-{}.json".format(s_id))

                with open(outFileJSON, 'r') as f:
                    data = json.load( f)

                if  data[s_id][u'TransMem'] == {}:
    #                 print 'No transmembrane resiudes^'
                    checked +=[s_id]
                else:
                    #the N/O membrane residues are not in the PDB file but OPM predicts we have a transmembrane protein
                    ends = []
                    for chain ,transmem_helix_list in data[s_id][u'TransMem'].items():
                        for transmem_helix in transmem_helix_list:
                            ends +=["{}_{}".format(chain, transmem_helix[0])]
                            ends += [ "{}_{}".format(chain, transmem_helix[-1])]

                    dfpoints = dfpdb[dfpdb.index.isin(ends)]
                    if len(dfpoints) >= 6:

                        xs =  dfpoints.x.tolist()
                        ys =  dfpoints.y.tolist()
                        zs =  dfpoints.z.tolist()

                        # do fit
                        tmp_A = []
                        tmp_b = []
                        for i in range(len(xs)):
                            tmp_A.append([xs[i], ys[i], 1])
                            tmp_b.append(zs[i])
                        b = np.matrix(tmp_b).T
                        A = np.matrix(tmp_A)

                        # Manual solution
                        fit = (A.T * A).I * A.T * b
                        errors = b - A * fit
                        residual = np.linalg.norm(errors)
                        #segregate the points
                        for index, row in dfpoints.iterrows():
                            point  = np.array([row.get('x'),row.get('y'), row.get('z')])
                            if np.dot((fit[0],fit[1],-1., fit[2]), (point[0],point[1],point[2],1)) >0:
                                N_coords.update({index : list(point) })
                            else:
                                O_coords.update({index : list(point)})


        if O_coords == {} or N_coords == {}: 
            if len(O_coords) + len(N_coords) == 0:
                label = 'Globular'
                print  (len(dfpdb) , '\t', s_id, '\t', len(N_coords),'\t',  len(O_coords) ,'\t', label    )    

            else:
                label = 'Associated'
                print ( len(dfpdb) , '\t', s_id, '\t', len(N_coords),'\t',  len(O_coords) ,'\t', label  )      

            for index in (dfpdb.index.tolist()):
                row = dfpdb.loc[index]
                closest_dict.update({index: label})
                embedded_dict.update({index: False})


        else:
            label = 'Membrane'
            print  (len(dfpdb) , '\t', s_id, '\t', len(N_coords),'\t',  len(O_coords) ,'\t', label   )     

            N_vector, O_vector = get_vectors_ON(N_coords, O_coords)  
#             print  (len(dfpdb) , '\t', s_id, '\t', N_vector,'\t', O_vector ,'\t', label   )     

            for index in (dfpdb.index.tolist()):
                row = dfpdb.loc[index]
                x = row.get('x')
                y = row.get('y')
                z = row.get('z')

                distances = []
        #             if N_vector != []:
                N_dist = protein_geometry.find_distance_q_to_n_plane(q = [x,y,z],n_plane=N_vector)
                distances += [N_dist]
        #             if O_vector != []:
                O_dist = protein_geometry.find_distance_q_to_n_plane(q = [x,y,z],n_plane=O_vector)
                distances += [O_dist]

        #             if len(distances) > 1:
                N_proj_res, N_unit_vector  = protein_geometry.project_q_into_plane(q = [x,y,z],n_plane=N_vector)
                NO_mem_dist = protein_geometry.find_distance_q_to_n_plane(q = N_proj_res,n_plane=O_vector)

                closest_leaf, furthest_leaf = get_closest_leaf(N_dist, O_dist)

                distances.sort()

                if distances[-1] > NO_mem_dist:
                    embedded = False
                else:
                    embedded = True
                #     print N_dist, O_dist, NO_mem_dist, embedded

                closest_dict.update({index: closest_leaf})
                embedded_dict.update({index: embedded})

        dfpdb['closest_leaf' ] = pd.Series(closest_dict)
        dfpdb['embedded' ] = pd.Series(embedded_dict)
        dfpdb.to_csv(op.join(OPMcsvFolder, 'OPMcsv-{}.csv'.format(s_id)))

        if len(N_coords) >=3 and len(O_coords) >=3:
            N_vector, O_vector = get_vectors_ON(N_coords, O_coords)
            planes_dict.update({s_id : {'N_vector' :list(N_vector) , 'O_vector' : list(O_vector)}})
            previousData_planes_dict.update({s_id : {'N_vector' :list(N_vector) , 'O_vector' : list(O_vector)}})
            
            
            
        checked +=[s_id]
        
        
        with open(json_planes_path ,'w') as f:
            json.dump(previousData_planes_dict, f)
    log.info("Saving membrane planes dictionary...\n\t{}".format(json_planes_path))
           
    return planes_dict, checked


def run_005D1_membraneQCQA(planes_dict, query = 'OPM'):
    dfangles = setUpDF_andGetAngleOPM(leaf_vectors_all=planes_dict,
                                                query = query)

    dfangles =  get_thickness_areaOPM(dfangles, query = query)
    return dfangles

    