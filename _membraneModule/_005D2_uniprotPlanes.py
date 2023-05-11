# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)

from .membraneMetrics import *
from .uniprot_TM import *


def run_005D2_uniprotResidues(potential_membrane_structures,
                            potential_membrane_complex,
                            dfalldata,
                            
                       ):

    checked = []
    cplx_residues2 = {}
    dfmembrane_residues_uniprot = pd.DataFrame()
    dfmem_index = 0
    
    for cplx in tqdm(potential_membrane_complex):
        if cplx in checked:
            continue
    #     if cplx != 'ABC-15-CPLX':
    #         continue
    #     print '\n', cplx
        dfcplx = dfalldata[dfalldata.cplx == cplx]

        structure_chain_uniprot_mapping  = []
        for gene in dfcplx.gene.unique():
            dfcplx_gene = dfcplx[dfcplx.gene == gene]
            for index, row in dfcplx_gene.iterrows():
                mappingId =  "{}&{}&{}&{}".format(row.structureId, row.sfileChainId, row.geneSeqId,row.gene)
                if mappingId not in structure_chain_uniprot_mapping:
                    structure_chain_uniprot_mapping+=[mappingId]


        uniprot_structure_residues = {}
        for mappingId in structure_chain_uniprot_mapping:
            geneId = mappingId.split('&')[-1]
            uniprotId = mappingId.split('&')[-2]
            sfileChainId = mappingId.split('&')[-3]
            structureId = mappingId.split('&{}&{}'.format(sfileChainId,uniprotId))[0]



            
            uniprotTextFile = op.join(qspaceDirs['UniprotSeqsDir'], '{}.txt'.format(uniprotId))
            uniprotFastaFile = op.join(qspaceDirs['UniprotSeqsDir'], '{}.fasta'.format(uniprotId))
            
            if not op.exists(uniprotTextFile):
                print ('uniprot File DNE : {}\t{}'.format(geneId, op.basename(uniprotTextFile)))
                continue
                
#             record = SeqIO.read(uniprotFastaFile, "fasta")
#                 print (geneSeqId)
#             seq = SeqProp(id =uniprotId,seq = str(record.seq), sequence_path =uniprotFastaFile,metadata_path=uniprotTextFile)
            
#             print (seq.sequence_path)
#             print (seq.metadata_path)
#             print (seq)
#             seq = SeqProp(id=uniprotId, metadata_file=uniprotTextFile)
#             seq.uniprot = uniprotId
#             tm_residues = seq.extract_uniprot_residues_from_uniprot_tm_dataframe()
#             tm_residues = seq.extract_uniprot_residues_from_uniprot_tm_dataframe()
            (dfuniMEM, TM_features) = parse_uniprot_TM_features_from_txt(uniprotTextFile) #functions from Uniprot_TM
            tm_residues = extract_residues_from_uniprot_TMdataframe(dfuniMEM)#functions from Uniprot_TM
            
            for leaf, leaf_residues in tm_residues.items():
                if leaf_residues == []:
                    continue
                if structureId not in uniprot_structure_residues:
                    uniprot_structure_residues.update({structureId : {leaf : {sfileChainId: leaf_residues}}})
                elif leaf not in uniprot_structure_residues[structureId]:
                    uniprot_structure_residues[structureId].update({leaf: {sfileChainId: leaf_residues}})
                elif sfileChainId not in uniprot_structure_residues[structureId][leaf]:
                    uniprot_structure_residues[structureId][leaf].update({sfileChainId : leaf_residues})
                else:
                    uniprot_structure_residues[structureId][leaf][sfileChainId] += leaf_residues

                for residue in leaf_residues:

                    dfmembrane_residues_uniprot.loc[dfmem_index,'cplx'] = cplx
                    dfmembrane_residues_uniprot.loc[dfmem_index,'structureId'] = structureId
                    dfmembrane_residues_uniprot.loc[dfmem_index,'geneId'] = geneId
                    dfmembrane_residues_uniprot.loc[dfmem_index,'sfileChainId'] = sfileChainId
                    dfmembrane_residues_uniprot.loc[dfmem_index,'uniprotId'] = uniprotId
    #                 dfmembrane_residues_uniprot.loc[dfmem_index,'mappingId'] = mappingId
                    dfmembrane_residues_uniprot.loc[dfmem_index,'leafId'] = leaf
                    dfmembrane_residues_uniprot.loc[dfmem_index,'uniprotResidueId'] = residue
                    dfmem_index +=1

            if tm_residues == {'Unknown_Residues': []}:
                continue
    #         print geneId, uniprotId, chainId, structureId,tm_residues


        if uniprot_structure_residues != {}:
            cplx_residues2.update({cplx : uniprot_structure_residues})
        checked += [cplx]

    # these are the uniprot sequence residues !!!!
    # we will need the residue # on the structure
    # still need to align these residue numbers to the structure reside numbers
    missing_resnum  = 0 
    
    for cplx in tqdm(dfmembrane_residues_uniprot.cplx.unique()):
        dfm_c  = dfmembrane_residues_uniprot[dfmembrane_residues_uniprot.cplx==cplx]
        dfa_c = dfalldata[dfalldata.cplx == cplx]

        for structureId in dfm_c.structureId.unique():
            dfm_c_s = dfm_c[dfm_c.structureId==structureId]
            dfa_c_s = dfa_c[dfa_c.structureId==structureId]


            for sfileChainId in dfm_c_s.sfileChainId.unique():
                dfm_c_s_ch = dfm_c_s[dfm_c_s.sfileChainId==sfileChainId]
                dfa_c_s_ch = dfa_c_s[dfa_c_s.sfileChainId==sfileChainId]

                for index, row in dfm_c_s_ch.iterrows():
                    dfa_c_s_ch_rN = dfa_c_s_ch[dfa_c_s_ch.seqNum == row.uniprotResidueId]
                    dfa_c_s_ch_rN = dfa_c_s_ch_rN[dfa_c_s_ch_rN.structNum.isna()== False]

                    if len(dfa_c_s_ch_rN) > 1:
                        log.warn("1 residue represented twice in structure |{}-{}-{}-{}-{}-{}".format(index, row.cplx, row.structureId, row.sfileChainId, row.uniprotId, row.leafId))
                        continue 
#                         raise (AttributeError,'Cannot have 1 residue represented by the same structure twice')

                    if len(dfa_c_s_ch_rN) == 0 :
                        missing_resnum +=1
                        log.warn("Missing residue |{}-{}-{}-{}-{}-{}".format(index, row.cplx, row.structureId, row.sfileChainId, row.uniprotId, row.leafId))
                        continue

                    structureResidueId = int(dfa_c_s_ch_rN.structNum.values[0])
                    dfmembrane_residues_uniprot.loc[index,'structureResidueId'] = structureResidueId
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '005D2-MembraneResiduesUniProt.json')
    log.info("Saving membrane residues dataframe ...\n\t{}".format(outfile))
    dfmembrane_residues_uniprot.to_csv(outfile)
    return dfmembrane_residues_uniprot

def fix_for_ecoli_proteome(dfmembrane_residues_uniprot):
#     print 'Fixing specific Uniprot -membrane leaflet annotations in E.coli....'

    
    fix_list = []

    # [Cytoplasmic, Extracellular] --> [Cytoplasmic, Periplasmic_C] 
    for uniprotId in ['P0AAC6','P0AC23','P0AEJ0','P0AFH2','P0AFH6','P0AAE8','P0ABB8']:
        dfuni = dfmembrane_residues_uniprot[dfmembrane_residues_uniprot.uniprotId == uniprotId]
        for index, row in dfuni.iterrows():
            if row.leafId == 'Extracellular':
                data = "{}&{}&{}&{}&{}&{}".format(row.cplx,row.structureId, row.uniprotId, row.sfileChainId, row.uniprotResidueId, 'Periplasmic_C')
                if data not in fix_list:
                    fix_list +=  [data]

    # [Periplasmic, Unknown_Residues] --> [Periplasmic_E, Extracellular]      
    for uniprotId in ['P0A940']:
        dfuni = dfmembrane_residues_uniprot[dfmembrane_residues_uniprot.uniprotId == uniprotId]
        for index, row in dfuni.iterrows():
            if row.leafId == 'Unknown_Residues':
                data = "{}&{}&{}&{}&{}&{}".format(row.cplx,row.structureId, row.uniprotId, row.sfileChainId, row.uniprotResidueId, 'Extracellular')
                if data not in fix_list:
                    fix_list +=  [data]

            if row.leafId == 'Periplasmic':
                data = "{}&{}&{}&{}&{}&{}".format(row.cplx,row.structureId, row.uniprotId, row.sfileChainId, row.uniprotResidueId, 'Periplasmic_E')
                if data not in fix_list:
                    fix_list +=  [data]
                    
    # [TM_1, TM_2] --> [Cytoplasmic, Periplasmic_C]      
    for uniprotId in ['P10906','P75799']:
        dfuni = dfmembrane_residues_uniprot[dfmembrane_residues_uniprot.uniprotId == uniprotId]
        for index, row in dfuni.iterrows():
            if row.leafId == 'TM_1':
                data = "{}&{}&{}&{}&{}&{}".format(row.cplx,row.structureId, row.uniprotId, row.sfileChainId, row.uniprotResidueId, 'Cytoplasmic')
                if data not in fix_list:
                    fix_list +=  [data]

            if row.leafId == 'TM_2':
                data = "{}&{}&{}&{}&{}&{}".format(row.cplx,row.structureId, row.uniprotId, row.sfileChainId, row.uniprotResidueId, 'Periplasmic_C')
                if data not in fix_list:
                    fix_list +=  [data]

    # [Unknown_Residues] --> [ Periplasmic_C]      
    wrong_input = 'Unknown_Residues'
    correct_output  = 'Periplasmic_C'
    for uniprotId in ['P0AGM2']:
        dfuni = dfmembrane_residues_uniprot[dfmembrane_residues_uniprot.uniprotId == uniprotId]
        for index, row in dfuni.iterrows():
            if row.leafId == wrong_input:
                data = "{}&{}&{}&{}&{}&{}".format(row.cplx,row.structureId, row.uniprotId, row.sfileChainId, row.uniprotResidueId, correct_output)
                if data not in fix_list:
                    fix_list +=  [data]
                
    fix_dict = {}
    for data in fix_list:
        leafId = data.split('&')[-1]
        uniprotResidueId = float(data.split('&')[-2])
        sfileChainId = data.split('&')[-3]
        uniprotId = data.split('&')[-4]
        structureId = data.split('&')[-5]
        cplxId = data.split('&{}&{}&{}'.format(structureId,uniprotId,sfileChainId))[0]

        df = dfmembrane_residues_uniprot[dfmembrane_residues_uniprot.cplx == cplxId]
        df = df[df.structureId == structureId]
        df= df[df.uniprotId == uniprotId]
        df = df[df.sfileChainId == sfileChainId]
        df = df[df.uniprotResidueId == uniprotResidueId]
        if len(df) != 1:
            log.warn("{}-{}-{}-{}-{}-{}".format(leafId, uniprotResidueId,sfileChainId,uniprotId,structureId,cplxId))
            continue
        index = df.first_valid_index()
        fix_dict.update({index : leafId})
    for index, leafId in fix_dict.items():
        dfmembrane_residues_uniprot.loc[index,'leafId'] = leafId
        
        
    #fix specific E.coli Residues
    infile = op.join(qspaceDirs['root_dir'], '_membraneModule/manualFix_UniProtresidues.json')
    if op.exists(infile):
        with open(infile, 'r') as f:
            specific_residues = json.load(f)
        print ('Fixing specific Uniprot residues in E.coli....\n\t> {}'.format(infile))

        for uniprotId , leafDict in specific_residues.items():
            for leafId, residueList in leafDict.items():
                for uniprotResidueId in residueList:
                    dfs = dfmembrane_residues_uniprot[dfmembrane_residues_uniprot.uniprotId ==uniprotId]
                    dfs = dfs[dfs.uniprotResidueId == uniprotResidueId]
                    if len(dfs) > 0:
                        print ('Fixed :\t' ,len(dfs), leafId, uniprotId , uniprotResidueId)
                    for index, row in dfs.iterrows():
                        dfmembrane_residues_uniprot.loc[index,"leafId"] = leafId
    else:
        print ('no specific residues')
    return dfmembrane_residues_uniprot


def run_005D2_uniprotPlanes(dfmembrane_residues_uniprot,
                           sfileDict):
    leaf_vectors_all = {}
    leaf_coords_all = {}
    
#     print "Getting Uniprot membrane data ...\n------------------------------------------------"

    for cplx in tqdm(dfmembrane_residues_uniprot.cplx.unique()):
        dfc = dfmembrane_residues_uniprot[dfmembrane_residues_uniprot.cplx == cplx]
        dfc=dfc[dfc.structureResidueId.isna() == False]
        log.info(cplx)
    #     if cplx != 'CPLX-165':
    #         continue
        for structureId in dfc.structureId.unique():

            log.info(structureId)
            dfs = dfc[dfc.structureId == structureId]
            sfile = sfileDict[structureId]
            log.info(sfile)

            for leafId in dfs.leafId.unique():
                if leafId =='Unknown_Residues':
                    continue
                log.info(leafId)

                dfs_leaf = dfs[dfs.leafId == leafId]
                dfs_leaf = dfs_leaf[dfs_leaf.structureResidueId.isna() == False]
                leaflet_residues = []
                
                for index, row in dfs_leaf.iterrows():
                    leaflet_residues += ['{}_{}'.format(row.sfileChainId, int(row.structureResidueId))]
    #                 print int(row.structureResidueId)
                
                leaf_coords = get_coords_of_membrane_crossing_residues(sfile, leaflet_residues)
                if leaf_coords != {}:

                    N_vector, N_residues = protein_geometry.lstsq_determine_membrane_surface_plane(pdb_coords_per_gene_leaf=leaf_coords)
                    if len(N_vector) == 0:
                        continue
                    if cplx not in leaf_vectors_all:
                        leaf_vectors_all.update({cplx : {structureId:{leafId : list(N_vector)}}})
                        leaf_coords_all.update({cplx : {structureId:{leafId : leaf_coords}}})
                    elif structureId not in leaf_vectors_all[cplx]:
                        leaf_vectors_all[cplx].update({structureId:{leafId : list(N_vector)}})
                        leaf_coords_all[cplx].update({structureId:{leafId : leaf_coords}})
                    elif leafId not in  leaf_vectors_all[cplx][structureId]:
                        leaf_vectors_all[cplx][structureId].update({leafId : list(N_vector)})
                        leaf_coords_all[cplx][structureId].update({leafId : leaf_coords})

        
        
    leaf_vectors_all_json = copy.deepcopy(leaf_vectors_all)
    for c, cinfo in leaf_vectors_all.items():
        for s, leafinfo in cinfo.items():
            for leaf, leafplane in leafinfo.items():
                new_leafplane = []

                for value in leafplane:
                    new_leafplane += [float(value)]
    #             print new_leafplane 

                leaf_vectors_all_json[c][s][leaf] = new_leafplane
        
    json_planes_path = op.join(qspaceDirs['DataOutput_dir'], '005D2-MembranePlanesUniProt.json')
    with open(json_planes_path ,'w') as f:
        json.dump(leaf_vectors_all_json, f)
    log.info("Saving membrane planes ...\n\t{}".format(json_planes_path))

    leaf_vectors_all = copy.deepcopy(leaf_vectors_all_json)
    return leaf_vectors_all


def run_005D2_membraneQCQA(planes_dict,dfmembrane_residues_uniprot, sfileDict, query = 'Uniprot', csvFolder =  qspaceDirs['UniprotCsvDir']):
#     print "Setting up membrane QCQA and calculating angle between planes ...\n------------------------------------------------"
    dfangles = setUpDF_andGetAngleUniProt(leaf_vectors_all=planes_dict,
                                                query = query)
#     print "Area and thickness membrane QCQA ...\n------------------------------------------------"
    
    dfangles =  get_thickness_areaUniProt(dfangles=dfangles,
                                                       sfileDict = sfileDict,
                                                       query = query,
                                                       UniprotcsvFolder  =csvFolder,
                                                      )
    for index, row in dfangles.iterrows():
        df = dfmembrane_residues_uniprot[dfmembrane_residues_uniprot.cplx == row.cplx]
        df= df[df.structureId == row.structureId]
        dfangles.loc[index,'UniprotIds'] = str( df.uniprotId.unique().tolist())
    for index, row in dfangles.iterrows():
        dfangles.loc[index, 'sfile'] = sfileDict[row.structureId]
        
    outfile = op.join(qspaceDirs['DataOutput_dir'], '005D2-{}_membrane_data.csv'.format(query))
    dfangles.to_csv(outfile)
    log.info("Saving membrane calculations...\n\t{}".format(outfile))
    return dfangles

