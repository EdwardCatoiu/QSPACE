# -*- coding: utf-8 -*-
from .utils import *
from .uniprot_TM import *


from .membraneMetrics import *

log = logging.getLogger(__name__)

def get_bulb(index,csvFolder, query = 'Uniprot'):
    found = False
    for filename in os.listdir(csvFolder):
        if query + "csv" + '-'+ index + '.csv' in filename:
            found = True
            break
    if found:
#         print index
#         print filename
        dfpdb =pd.read_csv(op.join(csvFolder, filename),index_col = 0 )
        dfbulb = dfpdb[dfpdb.embedded == False]
        return dfpdb,dfbulb
def filter_angle(df, angle_cutoff):
    df = df[df.angle <= angle_cutoff]
    return df

def filter_thickness(df, membrane_thickness = [12,45]):
    df = df[df.membrane_thickness >= membrane_thickness[0]]
    df = df[df.membrane_thickness <= membrane_thickness[1]]
    
    return df

def filter_area(df, area_upper_lim = 100000 ):
    embedded_areas = {}
    for index, row in df.iterrows():
        area_total = 0 
        is_membrane_embedded = True
        for leaf, area in ast.literal_eval(row.embedded_area).items():
            if area ==0:
#                 print index, row.structureId
                is_membrane_embedded = False

            area_total += area

    #     print "{:.0f}".format(area_total), '\t',index, row.structureId
        if is_membrane_embedded:
            embedded_areas.update({index : area_total/2.})

    df['avg_area_per_leaflet'] = pd.Series(embedded_areas)
    
    df = df[df.avg_area_per_leaflet.isna() == False]
    df = df[df.avg_area_per_leaflet <= area_upper_lim]
    
    return df

def rename_chains_OPM(index,dfpdb_opm,structureOPMChains):
    #this block is to modify in opm chains to reflect those that were changed before making the PDBs in 005B that were
    #then sent to OPM servers for calculations.
    if index in structureOPMChains:
        orig_chain_opm_chain = structureOPMChains[index]
        if type(orig_chain_opm_chain) == str:
            orig_chain_opm_chain = ast.literal_eval(orig_chain_opm_chain)
        for index2, row2 in dfpdb_opm.iterrows():
            opm_chain, opm_residue = row2.chain_residue.split('_')

            for orig_chain, opm_mod_chain in orig_chain_opm_chain.items():
                if opm_chain == opm_mod_chain:
                    #this is where we swap back the original chain into the chain_residue which mirrors the index
                    # for df opm pdb
                    new_chain_residue = "{}_{}".format(orig_chain, opm_residue)
                    dfpdb_opm.loc[index2,'chain_residue'] = new_chain_residue
#                         print index, index2, new_chain_residue
    return dfpdb_opm

def get_leaf_translation(df,uniprot_leaves,opm_leaves):
        leaf_translation = {}
        maximum_residues = df.max().max()
        while maximum_residues is not np.nan:

            maximum_residues = df.max().max()

            for index, row in df.iterrows():
                for column in df.columns:
                    if row.get(column) == maximum_residues:
#                         print index, column
                        df = df.drop(index)
                        df = df.drop(column, axis = 1)
                        leaf_translation.update({index : column})
        
        for opm_leaf, uniprot_leaf in leaf_translation.items():
            if opm_leaf + '_vector' in opm_leaves:
                opm_leaves.remove(opm_leaf + '_vector')
            if uniprot_leaf  in uniprot_leaves:
                uniprot_leaves.remove(uniprot_leaf)

        if len(uniprot_leaves) == len(opm_leaves) == 1:
            leaf_translation.update({opm_leaves[0].split('_vector')[0] : uniprot_leaves[0]})
        elif len(uniprot_leaves) > 1 or len(opm_leaves) > 1:
            raise ValueError
        
        return leaf_translation
    
##########################################################################
def run_005E_combineMembraneData(membrane_calculations):
    for source, df in membrane_calculations.items():
        dfout = filter_angle(df, angle_cutoff= 35)
        dfout = filter_thickness(dfout, membrane_thickness = [12,45])
        dfout = filter_area(dfout, area_upper_lim = 100000 )
#         print len(df), len(dfout), source

        membrane_calculations.update({source : dfout})
        
        
        
       
        
        
    structures = []
    for source, df in membrane_calculations.items():
        structures += df.structureId.tolist()
    structures = list(set(structures))
    
    ####################################################################
    dfsource = pd.DataFrame()
    for structureId in structures:
        for source, df in membrane_calculations.items():
            if structureId in df.structureId.values:
                leaves = df[df.structureId == structureId].leaf1.tolist() + df[df.structureId == structureId].leaf2.tolist()  
                leaves = list(set(leaves))

                dfsource.loc[structureId, source] = str(leaves)

    for index, row in dfsource.iterrows():
        if 'OPM' in row.dropna().keys():
            dfsource.loc[index,'final_source'] = 'OPM'
        elif 'Uniprot' in row.dropna().keys():
            dfsource.loc[index,'final_source'] = 'Uniprot'
        elif 'TMHMM' in row.dropna().keys():
            dfsource.loc[index,'final_source'] = 'TMHMM'

        if "Uniprot" in row.dropna().keys():
            if ast.literal_eval(row.get('Uniprot')) != ['TM_1','TM_2']:
                dfsource.loc[index,'orient_Uniprot'] = True
    
    
    uniprotIds_dict ={}
    dfuni  = membrane_calculations["Uniprot"]
    for structureId in dfuni.structureId.unique():
        dfs = dfuni[dfuni.structureId == structureId]
        uniprotIds = []
        for index, row in dfs.iterrows():
            uniprotIds += ast.literal_eval(row.UniprotIds)
    #     print structureId, uniprotIds
        uniprotIds_dict.update({structureId : uniprotIds})
    dfsource['UniprotIds'] = pd.Series(uniprotIds_dict)

    for structureId in ['5nik-assembly1','5v5s-assembly1','5o66-assembly1','5ng5-assembly1']:
        if structureId in dfsource.index:
            if 'Uniprot' in dfsource.loc[structureId].dropna().keys():
                 dfsource.loc[structureId,'final_source']= 'Uniprot'
    return membrane_calculations,dfsource
    
    
def run_005E_UniOPMshared_bulbs(dfsource,structureOPMChains):
    ##################################################################  
    shared_dict  = {}
    print ("Finding shared residues between OPM and UniProt bulbs\n--------------------------------------")
    for index in tqdm(dfsource[dfsource.final_source =='OPM'].index.tolist()):
        row = dfsource.loc[index]
    #     break
        if 'Uniprot' in row.dropna().keys() and 'OPM' in row.dropna().keys():
            dfpdb_uni,dfbulb_uni = get_bulb(index,csvFolder = qspaceDirs['UniprotCsvDir'],query = 'Uniprot')
            dfpdb_opm,dfbulb_opm = get_bulb(index,csvFolder = qspaceDirs['opmCsvDir'],query = 'OPM')
            dfpdb_opm['chain_residue'] = dfpdb_opm.index

            dfpdb_opm = rename_chains_OPM(index,dfpdb_opm,structureOPMChains)
            dfbulb_opm = dfpdb_opm[dfpdb_opm.embedded == False]

            dfshared = pd.DataFrame()
    #         print index
            for leaf1 in dfbulb_opm.closest_leaf.unique():
                dfbulb_opm_leaf = dfbulb_opm[dfbulb_opm.closest_leaf == leaf1]

                for leaf2 in dfbulb_uni.closest_leaf.unique():
                    dfbulb_uni_leaf = dfbulb_uni[dfbulb_uni.closest_leaf == leaf2]

                    shared = set(dfbulb_opm_leaf.chain_residue.tolist()).intersection(dfbulb_uni_leaf.index.tolist())

                    dfshared.loc[leaf1,leaf2] = len(shared)
            if dfshared.max().max() == 0:
                log.info("No shared bulb residues in OPM/Uniprot {}".format(index))
                continue
            shared_dict.update({index : dfshared})  
     
    return shared_dict
    
    ##################################################################  
    
def run_005E_auto_orient(shared_dict,dfsource):    
    leaf_translation = {}
    e = 0 
    for index, df in tqdm(shared_dict.items()):
#         print (len(df), len(df.columns), index)
    #     print index
        try:
            uniprot_leaves = ast.literal_eval(dfsource.loc[index,'Uniprot'])
            if set(uniprot_leaves) == set(['TM_1','TM_2']):
                continue
            uniprot_leaves = list(set(uniprot_leaves) - set(['TM_1','TM_2']))
            opm_leaves =  ast.literal_eval(dfsource.loc[index,'OPM'])
            leaf_translation.update({index : get_leaf_translation(df,uniprot_leaves,opm_leaves)}) 
    #     print index, leaf_translation
        except ValueError:
            log.warn("Error : {}".format(index))
            e +=1
            continue

    for index in ['5nik-assembly1','5v5s-assembly1','TORS-CPLX_AlphaMulti','5o66-assembly1','5ng5-assembly1',]:
        if index not in dfsource.index:
            continue
        row = dfsource.loc[index]
        leaves  = row.Uniprot
        if type(leaves) == str:
            leaves = ast.literal_eval(leaves)
#         print index, leaves
        leafdict=  {}
        for leaf in leaves:
            leafdict.update({leaf:leaf})
        leaf_translation.update({index : leafdict})
    dfsource['Uniprot_bulb_orient'] = pd.Series(leaf_translation)
    
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '005E-QCQA_dfsource.csv')
    log.info("Saving membrane calculations for QCQA'd proteins...\n\t{}".format(outfile))
    dfsource.to_csv(outfile)
    
    dfmissing = dfsource[dfsource.Uniprot_bulb_orient.isna()]
    dfmissing = dfmissing[dfmissing.OPM.isna() == False]
    
    return dfsource, dfmissing


def get_topoDF(geneId,uniprotId):
    uniprotTextFile = op.join(qspaceDirs['UniprotSeqsDir'], '{}.txt'.format(uniprotId))
    (dfuniMEM, TM_features) = parse_uniprot_TM_features_from_txt(uniprotTextFile) #functions from Uniprot_TM
#     tm_residues = extract_residues_from_uniprot_TMdataframe(dfuniMEM)#functions from Uniprot_TM
    
#     seq = SeqProp(ident=uniprotId, metadata_file=uniprotTextFile)
#     seq.uniprot = uniprotId
#     topological_domains  = seq.uniprot_tm_dataframe
    topological_domains  = dfuniMEM
    topological_domains=topological_domains[topological_domains.get('type') == 'TOPO_DOM']
    if set(list(topological_domains.location.unique()))== set(['Cytoplasmic','Periplasmic']):
        topological_domains = topological_domains.replace('Periplasmic','Periplasmic_C')

    elif set(list(topological_domains.location.unique()))== set(['Extracellular','Periplasmic']):
        topological_domains = topological_domains.replace('Periplasmic','Periplasmic_E')
          
    return topological_domains
 
def get_topoRES(topological_domains):

    topological_residues = {}
    for index, row in topological_domains.iterrows():
        for uniprotResNumber in range(int(row.seq_start),int(row.seq_end) + 1):
            topological_residues.update({uniprotResNumber :row.location})
    return topological_residues



def run_005E_gapFillMissingTopoDomains(dfsource,
                                       dfmissing,
                                       dfalldata,
                                       structureOPMChains,
                                       missingleaves={},
                                       leaf_translation_failedQCQA= {},
                                      ):
    
    shared_dict_failedQCQA= {}
    checked = []

    for index in tqdm(dfmissing.index.tolist()):
#         print (index)
        if index in checked:
            continue
        row = dfmissing.loc[index]
        dfpdb_opm,dfbulb_opm = get_bulb(index,csvFolder = qspaceDirs['opmCsvDir'],query = 'OPM')
        dfpdb_opm['chain_residue'] = dfpdb_opm.index
        dfpdb_opm = rename_chains_OPM(index,dfpdb_opm,structureOPMChains)
        dfpdb_opm = dfpdb_opm.set_index("chain_residue",drop = True)

        dfbulb_opm = dfpdb_opm[dfpdb_opm.embedded == False]

        opm_leaves = ast.literal_eval(row.OPM)

        dfs = dfalldata[dfalldata.structureId == index]
    #     shared_dict_uni = {}
        log.info(index)
        uniprotTopoLocation_dict = {}
        for uniprotId in dfs.geneSeqId.unique():
            log.info(uniprotId)
            dfu = dfs[dfs.geneSeqId == uniprotId]
            for chain in dfu.sfileChainId.unique():
                log.info(chain)

                dfchain = dfu[dfu.sfileChainId == chain]
                geneId = dfchain.gene.unique()[0]
                log.info(geneId)

    #             print index, uniprotId,chain, geneId

                ###################################################################################################
                try:
                    topological_domains = get_topoDF(geneId,uniprotId)
                    topological_domains = topological_domains[topological_domains.seq_start != '?']
                    topological_domains = topological_domains[topological_domains.seq_end != '?']
                    topological_residues= get_topoRES(topological_domains)

                
                except IOError:
                    continue
                topological_residues= get_topoRES(topological_domains)

                
                for uniprotResnum, uniprotTopoLocation in topological_residues.items():
                    try:
                        cr = dfchain[dfchain.seqNum == uniprotResnum].Chain_Residue.values[0]
                        uniprotTopoLocation_dict.update({cr : uniprotTopoLocation})
                    except IndexError:
                        log.warn("unable to map {} - {} - {}".format(uniprotResnum, uniprotTopoLocation,cr))
                        continue
        if uniprotTopoLocation_dict == {}:
            checked +=[index]
            continue
    #     print 'found^'
        dfpdb_opm['UniprotTopologicalDomain'] = pd.Series(uniprotTopoLocation_dict)
        df = dfpdb_opm[dfpdb_opm.UniprotTopologicalDomain.isna() == False]
        dfshared = pd.DataFrame()
        for opm_leaf in df.closest_leaf.unique():
            dfopmleaf = df[df.closest_leaf == opm_leaf]

            for uniprot_leaf in df.UniprotTopologicalDomain.unique():
                dfopm_uniprotleaf = df[df.UniprotTopologicalDomain == uniprot_leaf]
                shared = len(set(dfopm_uniprotleaf.index.tolist()).intersection(dfopmleaf.index.tolist()))
        #         print shared, opm_leaf, uniprot_leaf
                dfshared.loc[opm_leaf, uniprot_leaf] = shared

    #             if opm_leaf not in shared_dict_uni:
    #                 shared_dict_uni.update({opm_leaf : {uniprot_leaf : shared}})
    #             elif uniprot_leaf not in shared_dict_uni[opm_leaf]:
    #                 shared_dict_uni[opm_leaf].update({uniprot_leaf : shared})
    #             else:
    #                 shared_dict_uni[opm_leaf][uniprot_leaf]+=shared

        uniprot_leaves =  df.UniprotTopologicalDomain.dropna().unique().tolist()
#         print ('\t', uniprot_leaves, opm_leaves)
#         print (opm_leaves)
        is_ok = True
        if len(uniprot_leaves) == 1:
            if 'Cytoplasmic' in uniprot_leaves:
                uniprot_leaves +=['Periplasmic_C']
            elif 'Periplasmic_C' in uniprot_leaves:
                uniprot_leaves +=['Cytoplasmic']
            elif 'Extracellular' in uniprot_leaves:
                uniprot_leaves +=['Periplasmic_E']

            elif index in missingleaves:
                uniprot_leaves += missingleaves[index]
            else:
                log.warn('Error w/ {}-{}-{}'.format(index, uniprot_leaves, opm_leaves))
                is_ok = False
                
        if not is_ok:
            continue
        if len(uniprot_leaves) ==0:
            log.warn('No Uniprot Leaves can be mapped to proteome: Check uniprot.txt, or QSPACE mapping (structAA1)\n{}-{}-{}'.format(index, uniprot_leaves, opm_leaves))

            continue
            
#         log.info('{}-{}-{}'.format(index, uniprot_leaves, opm_leaves))

        leaf_translation_failedQCQA.update({index: get_leaf_translation(dfshared,uniprot_leaves,opm_leaves )})
        shared_dict_failedQCQA.update({index : dfshared})
        checked +=[index]
    


#     leaf_translation_failedQCQA.update({
#     #     "CADB_ECOLI_model1_clean_residues_removed" :{'O': 'Periplasmic_C', 'N': 'Cytoplasmic'},
#                                        "AF-P02929-F1-model_v2" :{'O': 'Periplasmic_C', 'N': 'Cytoplasmic'},
#     #                                    "EMRB_ECOLI_model1_clean_residues_removed" :{'O': 'Cytoplasmic', 'N': 'Periplasmic_C'},
#                                        "AF-P45762-F1-model_v2" : {'O': 'Periplasmic_C', 'N': 'Cytoplasmic'},
#     #     'P0AAC6_14_217_4pgs' : {'O': 'Periplasmic_C', 'N': 'Cytoplasmic'},
#     #     'P0ABB8_43_898_4hqj': {'O': 'Periplasmic_C', 'N': 'Cytoplasmic'},
#         'AF-P0AAE8-F1-model_v2' : {'O': 'Periplasmic_C', 'N': 'Cytoplasmic'},

#                                        }
#                                       )  
    dfsource['Uniprot_bulb_orient_failedQCQA'] = pd.Series(leaf_translation_failedQCQA)   

    keep = []
    for index, row in dfsource.iterrows():
        if set(row.dropna().keys()).intersection(set(['Uniprot_bulb_orient','Uniprot_bulb_orient_failedQCQA'])) == set():
            keep += [index]

    outfile = op.join(qspaceDirs['DataOutput_dir'], '005E-needs_manual_orientation_230322.csv')
    dfsource[dfsource.index.isin(keep)].to_csv(outfile)
    log.info("Saving protein structures that need manual orientation ...\n\t{}".format(outfile))

    return dfsource,dfsource[dfsource.index.isin(keep)]


def run_005E_useManualOrient(dfsource,dfuni, dftmhmm,membrane_calculations, base = qspaceDirs['opmManualOrientDir']):
    location_dict = {}
    opmnotgood = []
    import os
    for root, dirs, files in os.walk(base, topdown=False):
        for name in files:
            fullPath = os.path.join(root, name)
            structureBase = op.basename(fullPath)
            structureId = "_".join(structureBase.split('_')[1:]).split('.png')[0]
            if structureId not in dfsource.index:
                continue

            if "O-Extracellular" in fullPath:
                location = {u'O': u'Extracellular', u'N': u'Periplasmic_E'}
            elif 'N-Extracellular' in fullPath:
                location = {u'N': u'Extracellular', u'O': u'Periplasmic_E'}
            elif 'O-Periplasmic_C' in fullPath:
                location = {u'N': u'Cytoplasmic', u'O': u'Periplasmic_C'}
            elif 'N-Periplasmic_C' in fullPath:
                location = {u'O': u'Cytoplasmic', u'N': u'Periplasmic_C'}
            elif 'IM/No_info/' in fullPath:
                location = {u'O': u'IM_1', u'N': u'IM_2'}
            elif '/OM/' in fullPath:
                location = {u'O': u'OM_1', u'N': u'OM_2'}
            elif '/Cytoplasmic/' in fullPath:
                location = "Cytoplasm"
            elif '/Periplasmic/' in fullPath:
                location = "Periplasm"
            elif '/Extracellular/' in fullPath:
                location = "Extracellular"
            else:
                print ('OPM is not good\n>{}'.format(op.basename(fullPath)))
                opmnotgood +=[structureId]
                continue

    #         print fullPath
            data = fullPath.split('OPMmanualPNGs')[1].split(structureBase)[0].split('/')
    #         print data
            while '' in data:
                data.remove('')

            for l in ['O-Periplasmic_C','N-Periplasmic_C','O-Extracellular','N-Extracellular',]:
                if l in data:
                    data.remove(l)

            data = "&".join(data)

            location_dict.update({structureId : location})

    dfsource['Manual_bulb_orient']  = pd.Series(location_dict)

####################################################
    t = []
    location_dict = {}
    final_leafInfo = {}
    indexlist1 = []
    indexlist2 = []
    for index, row in dfsource.iterrows():
        data = row[['Uniprot_bulb_orient','Uniprot_bulb_orient_failedQCQA','Manual_bulb_orient']]
        data=data.dropna()

        final_leaves = False

        if 'IM_1' in str(data.values):
            location = 'Inner_Membrane'
        elif "Periplasmic_C" in str(data.values):
            location = 'Inner_Membrane'    
        elif "Periplasmic_E" in str(data.values):
            location = 'Outer_Membrane'
        elif str(['Cytoplasm']) == str(data.values):
            location = 'Cytoplasm'
        elif str(['Periplasm']) == str(data.values):
            location = 'Periplasm'
        elif str(['Extracellular']) == str(data.values):
            location = 'Extracellular'
        elif 'OM_1' in str(data.values):
            location = 'Outer_Membrane'
        elif index in ['YJGL_ECOLI_model1_clean_residues_removed','2wcd-assembly1']:
            location ='Inner Membrane (No Area)'
        elif index in ['MDTQ_ECOLI_model1_clean_residues_removed']:
            location = 'Outer Membrane (No Area)' 
        elif index in ['P39285_831_1086_4hw9']:
            location ='Inner_Membrane'
            final_leaves = {"tmhmm_leaflet_O":'Periplasmic_C',"tmhmm_leaflet_I":'Cytoplasmic'}
        elif index in ['P76045_22_301_5mwv']:
            location ='Outer_Membrane'
            final_leaves = {"TM_1":'Extracellular',"TM_2":'Periplasmic_E'}
        elif index in ['P39367_1_248_1nkv']:
            location ='Cytoplasm'
        elif index =='AF-P02929-F1-model_v2':
            final_leaves = {'O': 'Periplasmic_E', 'N': 'Extracellular'}
            location ='Outer Membrane'
        elif index in ['P77747_23_377_6rck']:
            final_leaves = {'N': 'Periplasmic_E', 'O': 'Extracellular'}
            location ='Outer Membrane'
        

        else:
            t += [index]
            log.warn("{} needs to by oriented manually, download the sfile and put it in the appropriate folder".format(index))
#             print (index, data.values,'location missing',final_leaves)
            continue


        if not final_leaves:

            if len(data.values) == 0:
                final_leaves = {}
    #             print '0\t', index, data.values
                indexlist1 +=[index]

            elif len(data.values) != 1:
                final_leaves = {}
    #             print '1\t', index, data.values
                indexlist2 +=[index]

            elif type(data.values[0]) == str:
                if data.values[0] in ['Periplasm','Cytoplasm','Extracellular']:
    #                 print '2', index, data.values
                    final_leaves = {}
            else:

                final_leaves = data.values[0]
                if len(final_leaves) > 3:
                    location = 'TransMembrane'
    #                 print '3', index, data.values

        final_leafInfo.update({index : final_leaves})
        location_dict.update({index : location})
    
    for structureId in ['7z0s-assembly1','7z7v-assembly1']:
        if structureId not in dfsource.index:
            continue
        final_leafInfo.update({structureId : {'O': 'Periplasmic_C', 'N': 'Cytoplasmic'}})

    dfsource['location'] = pd.Series(location_dict)
    dfsource['finalLeafInfo'] = pd.Series(final_leafInfo)
    
    
    ##############################################
    df3 = dfsource[dfsource.Manual_bulb_orient.isna() == False]
    
    dfsource_final = pd.DataFrame()
    dfsource_final['location'] = pd.Series(dfsource['location'])
    dfsource_final['finalLeafInfo'] = pd.Series(dfsource['finalLeafInfo'])
    dfsource_final['finalSource'] = pd.Series(dfsource['final_source'])


    for index, row in dfsource_final.iterrows():
        #is there uniprot data?
        if index in dfuni.structureId.unique():
            dfsource_final.loc[index,'Uniprot'] = 'Exists'
        if index in membrane_calculations['Uniprot'].structureId.unique():
            dfsource_final.loc[index,'Uniprot'] = 'Pass'

        #is there TMHMM data?
        if index in dftmhmm.structureId.unique():
            dfsource_final.loc[index,'TMHMM'] = 'Exists'
        if index in membrane_calculations['TMHMM'].structureId.unique():
            dfsource_final.loc[index,'TMHMM'] = 'Pass'

        #How were the leaves confirmed?
        if index in df3.index.unique():
            dfsource_final.loc[index,'orientBy'] = 'Manually' 
        if index in dfsource[dfsource.Uniprot_bulb_orient_failedQCQA.isna() == False].index.tolist():
            dfsource_final.loc[index,'orientBy'] = 'Uniprot_Failed'
        if index in dfsource[dfsource.Uniprot_bulb_orient.isna() == False].index.tolist():
            dfsource_final.loc[index,'orientBy'] = 'Uniprot_Passed'


    outfile = op.join(qspaceDirs['DataOutput_dir'], '005E-membrane_leaves_and_sources.csv')
    dfsource_final.to_csv(outfile)
    log.info("Saving membrane leaf orienation and sources of membrane information...\n\t{}".format(outfile))
    
    return dfsource,dfsource_final


def run_005E_mapOrientationToProteome(dforientation,
                                      dfalldata,
                                      structureOPMChains,
                                      OPMcsvFolder = qspaceDirs['opmCsvDir'],
                                      UniprotcsvFolder = qspaceDirs['UniprotCsvDir'],
                                      TMHMMcsvFolder = qspaceDirs['TMHMMcsvDir'],
                                     ):
    ############## change the leafs in the OPM/Uniprot/TMHMM .csv files #################################

    for index, row in tqdm(dforientation.iterrows()):
        if row.location not in ['Inner_Membrane','Outer_Membrane','TransMembrane','Extracellular']:
            continue

        row_finalLeafInfo = row.finalLeafInfo
        if type(row_finalLeafInfo) == str:
            row_finalLeafInfo = ast.literal_eval(row_finalLeafInfo) 

        if row.finalSource == 'OPM':
    #         continue
            csvfile = op.join(OPMcsvFolder, 'OPMcsv-{}.csv'.format(index))
            dfpdb = pd.read_csv(csvfile, index_col=0)
            dfpdb_copy = copy.deepcopy(dfpdb)



            for original_leaf, final_leaf in row_finalLeafInfo.items():
                dfpdb_copy = dfpdb_copy.replace(original_leaf,final_leaf)

            dfpdb['closest_leaf_oriented'] = pd.Series(dfpdb_copy.closest_leaf)
            dfpdb.to_csv(csvfile)
        elif row.finalSource =='Uniprot':
    #         break
            csvfile = op.join(UniprotcsvFolder, 'Uniprotcsv-{}.csv'.format(index))
            dfpdb = pd.read_csv(csvfile, index_col=0)
            dfpdb_copy = copy.deepcopy(dfpdb)

            for original_leaf, final_leaf in row_finalLeafInfo.items():
                dfpdb_copy = dfpdb_copy.replace(original_leaf,final_leaf)

            dfpdb['closest_leaf_oriented'] = pd.Series(dfpdb_copy.closest_leaf)
            dfpdb.to_csv(csvfile)

        elif row.finalSource == 'TMHMM':
            csvfile = op.join(TMHMMcsvFolder, 'TMHMMcsv-{}.csv'.format(index))
            dfpdb = pd.read_csv(csvfile, index_col=0)
            dfpdb_copy = copy.deepcopy(dfpdb)

            for original_leaf, final_leaf in row_finalLeafInfo.items():
                dfpdb_copy = dfpdb_copy.replace(original_leaf,final_leaf)

            dfpdb['closest_leaf_oriented'] = pd.Series(dfpdb_copy.closest_leaf)
            dfpdb.to_csv(csvfile)
        
    ############## get global orientation dict #################################
    global_orientation_dict = {}
    global_embedded_dict = {}
    
    for index, row in tqdm(dforientation.iterrows()):
        if row.location not in ['Inner_Membrane','Outer_Membrane','TransMembrane']:
            continue

        if row.finalSource == 'OPM':
            csvfile = op.join(OPMcsvFolder, 'OPMcsv-{}.csv'.format(index))
            dfpdb = pd.read_csv(csvfile, index_col=0)
            dfpdb['chain_residue'] = dfpdb.index
            dfpdb = rename_chains_OPM(index,dfpdb,structureOPMChains)
            dfpdb = dfpdb.set_index("chain_residue",drop = True)

        elif row.finalSource =='Uniprot':
            csvfile = op.join(UniprotcsvFolder, 'Uniprotcsv-{}.csv'.format(index))
            dfpdb = pd.read_csv(csvfile, index_col=0)

        elif row.finalSource == 'TMHMM':
            csvfile = op.join(TMHMMcsvFolder, 'TMHMMcsv-{}.csv'.format(index))
            dfpdb = pd.read_csv(csvfile, index_col=0)
#         print csvfile
        structureOrientation = dfpdb['closest_leaf_oriented'].to_dict()
        structureEmbeddedness = dfpdb['embedded'].to_dict()

        dfs = dfalldata[dfalldata.structureId == index]
        globalCR = dfs[dfs.sfileChain_Residue.isin(structureOrientation)].sfileChain_Residue.to_dict()

        for globIndex, Chain_Residue in globalCR.items():
            global_orientation_dict.update({globIndex : structureOrientation[Chain_Residue]})
            global_embedded_dict.update({globIndex : structureEmbeddedness[Chain_Residue]})

    outfile = op.join(qspaceDirs['DataOutput_dir'],'005E2-StructuralProteome_membrane_orientation.json')
    with open(outfile, 'w') as f:
        json.dump(global_orientation_dict, f)
    log.info("Saving AA-level membrane orientation ...\n\t{}".format(outfile))

    outfile = op.join(qspaceDirs['DataOutput_dir'],'005E2-StructuralProteome_membrane_embeddedness.json')
    with open(outfile, 'w') as f:
        json.dump(global_embedded_dict, f)
    log.info("Saving AA-level membrane embeddedness ...\n\t{}".format(outfile))
    
    return global_orientation_dict,global_embedded_dict

