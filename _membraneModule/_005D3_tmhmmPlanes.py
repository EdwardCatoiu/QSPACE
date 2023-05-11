# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)


from .membraneMetrics import *

def run_005D3_write_tmhmm_file(potential_membrane_complex,dfalldata,force_rerun):
    
    
    tmhmm_file_to_send = op.join(qspaceDirs['DataOutput_dir'], '005D3-to_TMHMM_server.txt')
    if op.exists(tmhmm_file_to_send) and not force_rerun:
        print ("using previous file")
        print ("SEND FILE : {} to {}.\nMay need to break it down into smaller files for DTU server".format(tmhmm_file_to_send,"https://dtu.biolib.com/DeepTMHMM"))
        return
    
    c = 0 
#     print "Writing sequences to ....\n\t> {}".format(tmhmm_file_to_send)
    with open(tmhmm_file_to_send, 'w') as f:

        for cplx in tqdm(potential_membrane_complex):
            dfcplx = dfalldata[dfalldata.cplx == cplx]

            structure_chain_uniprot_mapping  = []
            for gene in dfcplx.gene.unique():
                dfcplx_gene = dfcplx[dfcplx.gene == gene]
                for index, row in dfcplx_gene.iterrows():
                    mappingId =  "{}&{}&{}&{}".format(row.structureId, row.sfileChainId, row.geneSeqId,row.gene)
                    if mappingId not in structure_chain_uniprot_mapping:
                        structure_chain_uniprot_mapping+=[mappingId]

            for mappingId in structure_chain_uniprot_mapping:
                geneId = mappingId.split('&')[-1]
                uniprotId = mappingId.split('&')[-2]
                chainId = mappingId.split('&')[-3]
                structureId = mappingId.split('&{}&{}'.format(chainId,uniprotId))[0]

                uniprotTextFile = op.join(qspaceDirs['UniprotSeqsDir'], '{}.txt'.format(uniprotId))
                faa_file = op.join(qspaceDirs['UniprotSeqsDir'], '{}.fasta'.format(uniprotId))

                if not op.exists(uniprotTextFile):
                    print ('uniprot File DNE : {}\t{}'.format(geneId, op.basename(uniprotTextFile)))
                    continue
                record = SeqIO.read(faa_file, "fasta")
                seq = SeqProp(id =uniprotId,seq = str(record.seq), sequence_path =faa_file,metadata_path=uniprotTextFile)
                seq.uniprot = uniprotId


                f.write(">{}\n".format(mappingId))
                f.write("{}\n".format(seq.seq_str))
                c+=1

    c
    print ("SEND FILE : {} to {}.\nMay need to break it down into smaller files for DTU server".format(tmhmm_file_to_send,"https://dtu.biolib.com/DeepTMHMM"))
    


def getmeta(label):
    geneId = label.split(' | ')[0].split('&')[-1]
    uniprotId = label.split('&')[-2]
    chainId = label.split('&')[-3]
    structureId = label.split('&{}&{}&{}'.format(chainId, uniprotId,geneId))[0]
    designation = label.split(' | ')[-1]

    return geneId,uniprotId,chainId,structureId,designation    

def label_TM_tmhmm_residue_numbers_and_leaflets(tmhmm_seq):
    
    '''
    Determines the residue numbers of the TM-helix residues that cross the membrane and labels them by leaflet.
    
    Returns:
        leaflet_dict: a dictionary with leaflet_variable : [residue list] where the variable is inside or outside
        TM_boundary dict: outputs a dictionar with : TM helix number : [TM helix residue start , TM helix residue end]
    
    Args:
    
        tmhmm_seq = g.protein.representative_sequence.seq_record.letter_annotations['TM-tmhmm']
    
    '''
    TM_number_dict = {}
    T_index = []
    T_residue = []

    residue_count = 1
    for residue_label in tmhmm_seq:
        if residue_label == 'M':
            T_residue.append(residue_count)
        
        residue_count = residue_count + 1
    TM_number_dict.update({'T_residue' :T_residue})
    
    #finding the TM boundaries
    T_residue_list = TM_number_dict['T_residue']
    
    count = 0
    max_count = len(T_residue_list)-1
    TM_helix_count = 0
    TM_boundary_dict= {}

    while count <= max_count:
        #first residue = TM start
        if count == 0:
            TM_start = T_residue_list[count]
            count = count + 1
            continue
        #Last residue = TM end
        elif count == max_count:
            TM_end = T_residue_list[count]
            TM_helix_count = TM_helix_count + 1
            TM_boundary_dict.update({'TM_helix_'+ str(TM_helix_count): [TM_start, TM_end]})
            break
        #middle residues need to be start or end
        elif T_residue_list[count] != T_residue_list[count + 1] - 1:
            TM_end = T_residue_list[count]
            TM_helix_count = TM_helix_count + 1
            TM_boundary_dict.update({'TM_helix_'+ str(TM_helix_count): [TM_start, TM_end]})
            #new TM_start
            TM_start = T_residue_list[count + 1]
        count = count + 1
    #assign leaflet to proper TM residues O or I
    leaflet_dict = {}
    for leaflet in ['O', 'I']:
        leaflet_list = []
        for TM_helix, TM_residues in TM_boundary_dict.items():
            for residue_num in TM_residues:
                tmhmm_seq_index = residue_num - 1
                previous_residue = tmhmm_seq_index - 1
                next_residue = tmhmm_seq_index + 1
                #identify if the previous or next residue closest to the TM helix start/end is the proper leaflet
                if tmhmm_seq[previous_residue] == leaflet or tmhmm_seq[next_residue] == leaflet:
                    leaflet_list.append(residue_num)
        leaflet_dict.update({'tmhmm_leaflet_' + leaflet : leaflet_list})

    return TM_boundary_dict, leaflet_dict


def run_005D3_tmhmmResidues(dfalldata, tmhmmFolder = qspaceDirs['tmhmmResultsDir']):
    dfDeepTMHMM_results = pd.DataFrame()
    index = 0 
    
    tmhmm_structureId = {}
    tmhmm_geneId = {}
    tmhmm_uniprotId = {}
    tmhmm_chainId = {}
    tmhmm_designation = {}
    tmhmm_tmhmmSeq = {}
    tmhmm_uniprotSeq = {}

    
    
    print ("Reading deepTMHMM results....\n-----------------------------------------")
    for filename in tqdm(os.listdir(tmhmmFolder)):
        
        filepath = op.join(tmhmmFolder, filename)
        print (filepath)
        f = open(filepath)
        text = f.read()
        info = text.split('>')

        for structureInfo in tqdm(info):
            if structureInfo.count('\n') == 0:
                continue
            label = structureInfo.split('\n')[0]
            uniprotSeq = structureInfo.split('\n')[1]
            tmhmmSeq = structureInfo.split('\n')[2]
            geneId,uniprotId,chainId,structureId,designation = getmeta(label)
            if geneId not in dfalldata.gene.unique():
                continue  #only use results from genes that are in our query list
            
#             dfDeepTMHMM_results.loc[index,'structureId'] = structureId
#             dfDeepTMHMM_results.loc[index,'geneId'] = geneId
#             dfDeepTMHMM_results.loc[index,'uniprotId'] = uniprotId
#             dfDeepTMHMM_results.loc[index,'chainId'] = chainId
#             dfDeepTMHMM_results.loc[index,'designation'] = designation
#             dfDeepTMHMM_results.loc[index,'tmhmmSeq'] = tmhmmSeq
#             dfDeepTMHMM_results.loc[index,'uniprotSeq'] = uniprotSeq
            
            tmhmm_structureId.update({index :structureId})
            tmhmm_geneId.update({index :geneId})
            tmhmm_uniprotId.update({index :uniprotId})
            tmhmm_chainId.update({index :chainId})
            tmhmm_designation.update({index :designation})
            tmhmm_tmhmmSeq.update({index :tmhmmSeq})
            tmhmm_uniprotSeq.update({index :uniprotSeq})
            
            index +=1
    
    dfDeepTMHMM_results['structureId'] = pd.Series(tmhmm_structureId)
    dfDeepTMHMM_results['geneId'] = pd.Series(tmhmm_geneId)
    dfDeepTMHMM_results['uniprotId'] = pd.Series(tmhmm_uniprotId)
    dfDeepTMHMM_results['chainId'] = pd.Series(tmhmm_chainId)
    dfDeepTMHMM_results['tmhmmSeq'] = pd.Series(tmhmm_tmhmmSeq)
    dfDeepTMHMM_results['uniprotSeq'] = pd.Series(tmhmm_uniprotSeq)
    dfDeepTMHMM_results['designation'] = pd.Series(tmhmm_designation)
    
    
    dfmembrane_index = 0
    print ("Adding Gene-Chain-Structure-Uniprot-Leaf information....\n-----------------------------------------")

    tmhmm_results_structureId={}
    tmhmm_results_geneId={}
    tmhmm_results_chainId={}
    tmhmm_results_uniprotId={}
    tmhmm_results_leafId={}
    tmhmm_results_uniprotResidueId={}
    
    
    for index in tqdm(dfDeepTMHMM_results.index.tolist()):
        row = dfDeepTMHMM_results.loc[index]
        (tm_helix, leafletDict )= label_TM_tmhmm_residue_numbers_and_leaflets(row.tmhmmSeq)
        for leaf, residues in leafletDict.items():
    #         print len(residues),
            for res in residues:
                tmhmm_results_structureId.update({dfmembrane_index: row.structureId})
                tmhmm_results_geneId.update({dfmembrane_index: row.geneId})
                tmhmm_results_chainId.update({dfmembrane_index:row.chainId})
                tmhmm_results_uniprotId.update({dfmembrane_index: row.uniprotId})
                tmhmm_results_leafId.update({dfmembrane_index:leaf})
                tmhmm_results_uniprotResidueId.update({dfmembrane_index: res})

                dfmembrane_index +=1
                
    dfmembrane_residues_tmhmm = pd.DataFrame()
    dfmembrane_residues_tmhmm['structureId'] = pd.Series(tmhmm_results_structureId)
    dfmembrane_residues_tmhmm['geneId'] = pd.Series(tmhmm_results_geneId)
    dfmembrane_residues_tmhmm['chainId'] = pd.Series(tmhmm_results_chainId)
    dfmembrane_residues_tmhmm['uniprotId'] = pd.Series(tmhmm_results_uniprotId)
    dfmembrane_residues_tmhmm['leafId'] = pd.Series(tmhmm_results_leafId)
    dfmembrane_residues_tmhmm['uniprotResidueId'] = pd.Series(tmhmm_results_uniprotResidueId)

    print ("Adding structure residues numbers....\n-----------------------------------------")

    checked = []      
    for index in tqdm(dfmembrane_residues_tmhmm.index.tolist()):
        if index in checked:
            continue
        row  = dfmembrane_residues_tmhmm.loc[index]
        df = dfalldata[dfalldata.structureId == row.structureId]
        df = df[df.sfileChainId == row.chainId]
        df = df[df.seqNum == row.uniprotResidueId  ]
        df = df[df.structNum.isna()== False]
        if len(df) > 1 and len(df.cplx.unique()) != len(df):
    #         raise AttributeError ('Cannot have 1 residue represented by the same structure twice')
            print ('Cannot have 1 residue represented by the same structure twice : ', index)
            continue

        if len(df) == 0 :
#             print index,
            continue
        structureResidueId = int(df.structNum.values[0])
        dfmembrane_residues_tmhmm.loc[index,'structureResidueId'] = structureResidueId
        checked+=[index]
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '005D3-MembraneResiduesTMHMM.csv')
#     print "Saving ... \n\t> {}".format(outfile)
    dfmembrane_residues_tmhmm.to_csv(outfile)
    return dfmembrane_residues_tmhmm

def run_005D3_tmhmmPlanes(dfmembrane_residues_tmhmm,
                          sfileDict,
                          dfalldata,
                         ):

    leaf_vectors_all = {}
    leaf_coords_all = {}



    for structureId in tqdm(dfmembrane_residues_tmhmm.structureId.unique()):
        if structureId not in dfalldata.structureId.unique():
            continue
        if structureId not in sfileDict:
            continue

        dfs = dfmembrane_residues_tmhmm[dfmembrane_residues_tmhmm.structureId == structureId]
        sfile = sfileDict[structureId]

        for leafId in dfs.leafId.unique():
            if leafId =='Unknown_Residues':
                continue
            dfs_leaf = dfs[dfs.leafId == leafId]
            dfs_leaf = dfs_leaf[dfs_leaf.structureResidueId.isna() == False]
            leaflet_residues = []
            for index, row in dfs_leaf.iterrows():
                leaflet_residues += ['{}_{}'.format(row.chainId, int(row.structureResidueId))]
    #                 print int(row.structureResidueId)
            leaflet_residues = list(set(leaflet_residues))
            leaf_coords = get_coords_of_membrane_crossing_residues(sfile, leaflet_residues)
            if leaf_coords != {}:

                N_vector, N_residues = protein_geometry.lstsq_determine_membrane_surface_plane(pdb_coords_per_gene_leaf=leaf_coords)
                if len(N_vector) == 0:
                    continue
                if structureId not in leaf_vectors_all:
                    leaf_vectors_all.update({structureId:{leafId : N_vector}})
                    leaf_coords_all.update({structureId:{leafId : leaf_coords}})
                elif leafId not in  leaf_vectors_all[structureId]:
                    leaf_vectors_all[structureId].update({leafId : N_vector})
                    leaf_coords_all[structureId].update({leafId : leaf_coords})
    
    
    leaf_vectors_all_json = copy.deepcopy(leaf_vectors_all)
    for s, leafinfo in leaf_vectors_all.items():
        for leaf, leafplane in leafinfo.items():
            new_leafplane = []

            for value in leafplane:
                new_leafplane += [float(value)]

            leaf_vectors_all_json[s][leaf] = new_leafplane

    json_planes_path = op.join(qspaceDirs['DataOutput_dir'], '005D3-MembranePlanesTMHMM.json')
    with open(json_planes_path ,'w') as f:
        json.dump(leaf_vectors_all_json, f)
#     print "Saving ... > {}\n".format(json_planes_path)

    leaf_vectors_all = copy.deepcopy(leaf_vectors_all_json)
    return leaf_vectors_all


def run_005D3_membraneQCQA(planes_dict, sfileDict, dfalldata, query = 'TMHMM', csvFolder =  qspaceDirs['TMHMMcsvDir']):
#     print "Setting up membrane QCQA and calculating angle between planes ...\n------------------------------------------------"
    dfangles = setUpDF_andGetAngleTMHMM(leaf_vectors_all=planes_dict,
                                                     query = query,
                                                     dfalldata = dfalldata, )
#     print "Area and thickness membrane QCQA ...\n------------------------------------------------"
   
    
    dfangles = get_thickness_areaTMHMM(dfangles=dfangles,
                                       sfileDict=sfileDict,
                                       query = 'TMHMM',
                                       TMHMMcsvFolder  = csvFolder,
                                      )
    
    return dfangles

                
