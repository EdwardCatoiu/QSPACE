from .utils import *


def setUpDF_andGetAngleOPM(leaf_vectors_all, query = 'OPM'):
    
    if query =='OPM':
        pairs = [['N_vector', 'O_vector']]
    
    c = 0 
    n = 0 
    dfangles  = pd.DataFrame()
    dfangles_index = 0 
    print ("Setting up membrane QCQA and calculating angle between planes ...\n------------------------------------------------")
    for structureId, leafInfo in leaf_vectors_all.items():
        for [leaf1,leaf2] in pairs:
            if leaf1 in leafInfo.keys() and leaf2 in leafInfo.keys():
                v1 = leafInfo[leaf1]
                v2 = leafInfo[leaf2] 

                u1 = protein_geometry.find_unit_norm(v1)
                u2 = protein_geometry.find_unit_norm(v2)

                angle = protein_geometry.find_angle_between_vector(u1,u2)
    #                 print structureId, angle, leaf1, leaf2
    #             dfangles.loc[dfangles_index, 'cplx' ] = cplx
                dfangles.loc[dfangles_index, 'structureId' ] = structureId
                dfangles.loc[dfangles_index, 'leaf1' ] = leaf1
                dfangles.loc[dfangles_index, 'leaf2' ] = leaf2
                dfangles.loc[dfangles_index, 'v1' ] = str(list(v1))
                dfangles.loc[dfangles_index, 'v2' ] = str(list(v2))
                dfangles.loc[dfangles_index, 'angle' ] = angle
                dfangles_index +=1
                
    return dfangles


def setUpDF_andGetAngleUniProt(leaf_vectors_all, query = 'Uniprot'):
    
    if query =='Uniprot':
        pairs = [['TM_1', 'TM_2'],['Periplasmic_C', 'Cytoplasmic'],['Extracellular', 'Periplasmic_E']]

    dfangles  = pd.DataFrame()
    dfangles_index = 0 
    print ("Setting up membrane QCQA and calculating angle between planes ...\n------------------------------------------------")

    for cplx, structures in leaf_vectors_all.items():
        for structureId, leafInfo in structures.items():
            found = False
            for [leaf1,leaf2] in pairs:
                if leaf1 in leafInfo.keys() and leaf2 in leafInfo.keys():
                    found = True
                    v1 = leafInfo[leaf1]
                    v2 = leafInfo[leaf2] 

                    u1 = protein_geometry.find_unit_norm(v1)
                    u2 = protein_geometry.find_unit_norm(v2)

                    angle = protein_geometry.find_angle_between_vector(u1,u2)
    #                 print structureId, angle, leaf1, leaf2
                    dfangles.loc[dfangles_index, 'cplx' ] = cplx
                    dfangles.loc[dfangles_index, 'structureId' ] = structureId
                    dfangles.loc[dfangles_index, 'leaf1' ] = leaf1
                    dfangles.loc[dfangles_index, 'leaf2' ] = leaf2
                    dfangles.loc[dfangles_index, 'v1' ] = str(list(v1))
                    dfangles.loc[dfangles_index, 'v2' ] = str(list(v2))
                    dfangles.loc[dfangles_index, 'angle' ] = angle
                    dfangles_index +=1
    return dfangles
                    
def setUpDF_andGetAngleTMHMM(leaf_vectors_all, dfalldata, query = 'TMHMM'):
 
    pairs = [['tmhmm_leaflet_O', 'tmhmm_leaflet_I']]

    c = 0 
    n = 0 
    dfangles  = pd.DataFrame()
    dfangles_index = 0 
    print ("Setting up membrane QCQA and calculating angle between planes ...\n------------------------------------------------")

    for structureId, leafInfo in tqdm(leaf_vectors_all.items()):
        found = False
        cplx = dfalldata[dfalldata.structureId == structureId].cplx.values[0]
        for [leaf1,leaf2] in pairs:
            if leaf1 in leafInfo.keys() and leaf2 in leafInfo.keys():
                found = True
                v1 = leafInfo[leaf1]
                v2 = leafInfo[leaf2] 

                u1 = protein_geometry.find_unit_norm(v1)
                u2 = protein_geometry.find_unit_norm(v2)

                angle = protein_geometry.find_angle_between_vector(u1,u2)

                dfangles.loc[dfangles_index, 'cplx' ] = cplx
                dfangles.loc[dfangles_index, 'structureId' ] = structureId
                dfangles.loc[dfangles_index, 'leaf1' ] = leaf1
                dfangles.loc[dfangles_index, 'leaf2' ] = leaf2
                dfangles.loc[dfangles_index, 'v1' ] = str(list(v1))
                dfangles.loc[dfangles_index, 'v2' ] = str(list(v2))
                dfangles.loc[dfangles_index, 'angle' ] = angle
                dfangles_index +=1
    return dfangles
    
                
def get_thickness_areaUniProt(dfangles,sfileDict, query = 'Uniprot',UniprotcsvFolder  = qspaceDirs['UniprotCsvDir']):
    print ("Area and thickness membrane QCQA ...\n------------------------------------------------")
    for cplx in tqdm(dfangles.cplx.unique()):
        dfc = dfangles[dfangles.cplx == cplx]

        for structureId in dfc.structureId.unique():
            sfile = sfileDict[structureId]
            dfs = dfc[dfc.structureId == structureId]

            vectorsAll={}
            for index, row in dfs.iterrows():
                vector_leaf1 = ast.literal_eval(row.get('v1'))
                vector_leaf2 = ast.literal_eval(row.get('v2'))
                leaf1 = row.get('leaf1')
                leaf2 = row.get('leaf2')
                vectorsAll.update( {leaf1 : vector_leaf1, leaf2: vector_leaf2})

            csvOutFile = op.join(UniprotcsvFolder,'{}csv-{}.csv'.format(query,structureId))
#             print csvOutFile
            if not op.exists(csvOutFile):
                dfpdb = get_pdb_coordinates(structureId,sfile)
                dfpdb = protein_geometry.closest_leaf_embeddedness_NEW(dfpdb, vectorsAll)
                dfpdb.to_csv(csvOutFile)

            else:
                dfpdb = pd.read_csv(csvOutFile, index_col = 0)
                dfpdb = protein_geometry.closest_leaf_embeddedness_NEW(dfpdb, vectorsAll)
                dfpdb.to_csv(csvOutFile)
    #         print '\n'
    #         print vectorsAll
    #         print csvOutFile

            for index, row in dfs.iterrows():
    #             print index
                vectors = {}
                vector_leaf1 = ast.literal_eval(row.get('v1'))
                vector_leaf2 = ast.literal_eval(row.get('v2'))
                leaf1 = row.get('leaf1')
                leaf2 = row.get('leaf2')
                vectors.update( {leaf1 : vector_leaf1, leaf2: vector_leaf2})
    #             print '\n', vectors

                #non_embedded_area
                leaf_polygons = protein_geometry.get_projected_polygons_to_leaf(dfpdb, vectors, only_use_embedded = False)
                non_embedded_area = {}
                for leaf, polygon_info in leaf_polygons.items():
                    polygon_coords= list(polygon_info.values())
                    if len(polygon_coords) < 3:
                        non_embedded_area.update({leaf : 0})
                        continue

                    polygon_vert_coords = protein_geometry.get_polygon_vertices(polygon_coords)
                    area = protein_geometry.find_area(np.array(polygon_vert_coords))
                    non_embedded_area.update({leaf : area})
    #             print  non_embedded_area

                #non_embedded_area
                leaf_polygons = protein_geometry.get_projected_polygons_to_leaf(dfpdb, vectors, only_use_embedded = True)
                embedded_area = {}
                for leaf, polygon_info in leaf_polygons.items():
                    polygon_coords= list(polygon_info.values())
                    if len(polygon_coords) < 3:
                        embedded_area.update({leaf : 0})
                        continue
                    polygon_vert_coords = protein_geometry.get_polygon_vertices(polygon_coords)
                    area = protein_geometry.find_area(np.array(polygon_vert_coords))
                    embedded_area.update({leaf : area})
    #             print embedded_area
                dfangles.loc[index, 'embedded_area'] = str(embedded_area)
                dfangles.loc[index, 'blub_area'] = str(non_embedded_area)
                dfangles.loc[index, 'membrane_thickness'] = dfpdb.membrane_thickness.mean()
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '005D2-{}_membrane_data.csv'.format(query))
    dfangles.to_csv(outfile)
    log.info("Saving membrane calculations...\n\t{}".format(outfile))
    return dfangles


def get_thickness_areaTMHMM(dfangles,sfileDict, query = 'TMHMM',TMHMMcsvFolder  = qspaceDirs['TMHMMcsvDir']):
    print ("Area and thickness membrane QCQA ...\n------------------------------------------------")

    for index in tqdm(dfangles.index.tolist()):
        row  = dfangles.loc[index]

        structureId = row.structureId
        sfile = sfileDict[structureId]


        vector_leaf1 = ast.literal_eval(row.get('v1'))
        vector_leaf2 = ast.literal_eval(row.get('v2'))
        leaf1 = row.get('leaf1')
        leaf2 = row.get('leaf2')

        vectors = {leaf1 : vector_leaf1, leaf2: vector_leaf2}


        tmhmm_outfile = op.join(TMHMMcsvFolder, '{}csv-{}.csv'.format(query, structureId))

        if not op.exists(tmhmm_outfile):
            dfpdb = get_pdb_coordinates(structureId,sfile)
            dfpdb = protein_geometry.closest_leaf_embeddedness_NEW(dfpdb, vectors)
            dfpdb.to_csv(tmhmm_outfile)

        else:
            dfpdb = pd.read_csv(tmhmm_outfile, index_col = 0)
            dfpdb = protein_geometry.closest_leaf_embeddedness_NEW(dfpdb, vectors)
            dfpdb.to_csv(tmhmm_outfile)

        #non_embedded_area
        leaf_polygons = protein_geometry.get_projected_polygons_to_leaf(dfpdb, vectors, only_use_embedded = False)
        non_embedded_area = {}
        for leaf, polygon_info in leaf_polygons.items():
            polygon_coords= list(polygon_info.values())
            if len(polygon_coords) < 3:
                non_embedded_area.update({leaf : 0})
                continue

            polygon_vert_coords = protein_geometry.get_polygon_vertices(polygon_coords)
            area = protein_geometry.find_area(np.array(polygon_vert_coords))
            non_embedded_area.update({leaf : area})

        #non_embedded_area
        leaf_polygons = protein_geometry.get_projected_polygons_to_leaf(dfpdb, vectors, only_use_embedded = True)
        embedded_area = {}
        for leaf, polygon_info in leaf_polygons.items():
            polygon_coords= list(polygon_info.values())
            if len(polygon_coords) < 3:
                embedded_area.update({leaf : 0})
                continue
            polygon_vert_coords = protein_geometry.get_polygon_vertices(polygon_coords)
            area = protein_geometry.find_area(np.array(polygon_vert_coords))
            embedded_area.update({leaf : area})

        dfangles.loc[index, 'embedded_area'] = str(embedded_area)
        dfangles.loc[index, 'blub_area'] = str(non_embedded_area)
        dfangles.loc[index, 'membrane_thickness'] = dfpdb.membrane_thickness.mean()
        dfangles.loc[index, 'sfile'] = sfile  
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '005D3-{}_membrane_data.csv'.format(query))
    dfangles.to_csv(outfile)
    log.info("Saving membrane calculations...\n\t{}".format(outfile))
    return dfangles


def get_thickness_areaOPM(dfangles, query = 'OPM'):
    print ("Area and thickness membrane QCQA ...\n------------------------------------------------")

    for index in tqdm(dfangles.index.tolist()):
        row  = dfangles.loc[index]
#         print  row.structureId

        structureId = row.structureId
    #     sfile = sfileDict[structureId]
        
        if query =='OPM':
            sfile = op.join(qspaceDirs['opmOutputStructuresDir'], "OPMstructure-{}out.pdb".format(structureId.replace('-','_')))
            csvFile = op.join(qspaceDirs['opmCsvDir'], 'OPMcsv-{}.csv'.format(structureId))


    #     sfile = '../../../../../mnt/wwn-0x5000c500b96b98e1-part1/Projects/E_coli/GEMPRO_220519/GEMPRO_220519/OPMresults/OPMstructures/OPMstructure-{}out.pdb'.format(structureId.replace('-','_'))

        vector_leaf1 = ast.literal_eval(row.get('v1'))
        vector_leaf2 = ast.literal_eval(row.get('v2'))
        leaf1 = row.get('leaf1')
        leaf2 = row.get('leaf2')

        vectors = {leaf1 : vector_leaf1, leaf2: vector_leaf2}


        if op.exists(csvFile):
            dfpdb = pd.read_csv(csvFile, index_col= 0)
        else:
            dfpdb = get_pdb_coordinates(structureId,sfile)


        dfpdb = protein_geometry.closest_leaf_embeddedness(dfpdb, leaf1,leaf2,vectors)

        #non_embedded_area
        leaf_polygons = protein_geometry.get_projected_polygons_to_leaf(dfpdb, vectors, only_use_embedded = False)
        non_embedded_area = {}
        for leaf, polygon_info in leaf_polygons.items():
            polygon_coords= list(polygon_info.values())
            if len(polygon_coords) < 3:
                non_embedded_area.update({leaf : 0})
                continue

            polygon_vert_coords = protein_geometry.get_polygon_vertices(polygon_coords)
            area = protein_geometry.find_area(np.array(polygon_vert_coords))
            non_embedded_area.update({leaf : area})

        #non_embedded_area
        leaf_polygons = protein_geometry.get_projected_polygons_to_leaf(dfpdb, vectors, only_use_embedded = True)
        embedded_area = {}
        for leaf, polygon_info in leaf_polygons.items():
            polygon_coords= list(polygon_info.values())
            if len(polygon_coords) < 3:
                embedded_area.update({leaf : 0})
                continue
            polygon_vert_coords = protein_geometry.get_polygon_vertices(polygon_coords)
            area = protein_geometry.find_area(np.array(polygon_vert_coords))
            embedded_area.update({leaf : area})

        dfangles.loc[index, 'embedded_area'] = str(embedded_area)
        dfangles.loc[index, 'blub_area'] = str(non_embedded_area)
        dfangles.loc[index, 'membrane_thickness'] = dfpdb.membrane_thickness.mean()
        dfangles.loc[index, 'sfile'] = sfile
    outfile = op.join(qspaceDirs['DataOutput_dir'], '005D1-OPM_membrane_data.csv')
    dfangles.to_csv(outfile)
    log.info("Saving membrane calculations...\n\t{}".format(outfile))
    return dfangles

    