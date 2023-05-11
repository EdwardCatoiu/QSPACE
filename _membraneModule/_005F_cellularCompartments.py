# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)


from .membraneMetrics import *



def run_005F_assignCompartmentsEcocycAndGO(dfalldata,
                                           global_embedded_dict,
                                           global_orientation_dict,
                                           dfrepseq,
                                           ecocyc_file = op.join(qspaceDirs['Input_dir'], '005A-EcocycSmartTable-All-genes-of-E.-coli-K-12-substr.-MG1655.txt'),        
                                          ):

    dfalldata['005E_Membrane_Embedded'] = pd.Series(global_embedded_dict)
    dfalldata['005E_Membrane_Orientation'] = pd.Series(global_orientation_dict)
    
    membraneStructures = dfalldata[dfalldata.get('005E_Membrane_Orientation').isna() == False].structureId.unique().tolist()
    non_membraneStructures = dfalldata[dfalldata.structureId.isin(membraneStructures) == False].structureId.unique().tolist()
    
    dfrepseq_uni =  dfrepseq[dfrepseq.UniProtId.isna() == False]
    
    dforientation = pd.DataFrame()
    genelist = {}
    uniprot_ids = {}
    for structureId in tqdm(non_membraneStructures):
        dfs =  dfalldata[dfalldata.structureId == structureId]
        genes = dfs.gene.unique().tolist()
        uniprotIds  =dfrepseq_uni[dfrepseq_uni.index.isin(genes)].UniProtId.tolist()
        genelist.update({structureId : genes})
        uniprot_ids.update({structureId : uniprotIds})
    dforientation['genelist'] = pd.Series(genelist)
    dforientation['uniprotIds']= pd.Series(uniprot_ids)
    
    ################### Gene Ontology ############################
    GO_meta = {}
    GO_meta_comp = {}
    for index, row in tqdm(dforientation.iterrows()):
        data = ''
        data_c = ''
        for uniprotId in row.uniprotIds:
            txtFile = op.join(qspaceDirs['UniprotSeqsDir'],'{}.txt'.format(uniprotId))
            if not op.exists(txtFile):
                print (txtFile)
                continue

            p, f, c= read_uniprot_text.find_GO(txtFile)
            data += str(p).lower() + str(f).lower() 
            data_c += str(c).lower()

        GO_meta.update({index : data})
        GO_meta_comp.update({index : data_c})
    
    dforientation['GO_metadata']= pd.Series(GO_meta)
    dforientation['GO_metadata_compartment']= pd.Series(GO_meta_comp)
    
    ################### ECOCYC ############################
    dfecocyc_loc = pd.read_csv(ecocyc_file, sep = '\t')
    dfecocyc_loc['gene'] = dfecocyc_loc['Accession-1']
    dfecocyc_loc = dfecocyc_loc[['Gene Name','gene','Locations']]
    dfecocyc_loc = dfecocyc_loc[dfecocyc_loc.gene.isna()==False]
    ECOCYC_loc = {}
    for index, row in dfecocyc_loc.iterrows():
        ECOCYC_loc.update({row.gene : row.Locations})
        
    ECO_meta = {}
    for index, row in tqdm(dforientation.iterrows()):
        data = ''
        for gene in row.genelist:
            if gene in ECOCYC_loc:
                data += str(ECOCYC_loc[gene])
        ECO_meta.update({index : data})

    dforientation['ECO_metadata_compartment']= pd.Series(ECO_meta)

    # determine metadata keywords for compartments
    for index, row in dforientation.iterrows():
        data = str(row.GO_metadata_compartment)+ str(row.ECO_metadata_compartment)
        found = False
        for kw in ['outer membrane','membrane-bounded', 'periplasmic space',
                   'cytosol','cytoplasm','membrane','periplasm',
                   'inner membrane','transmembrane','extracellular','cell wall']:
            if kw in data:
                dforientation.loc[index,kw]= True
                found = True
        if found:
            continue
#         print index, data
        
    # assign comparments using keywords    
    c = 1
    indexmissing = []
    for index, row in dforientation.iterrows():
        kwList = set(row.drop(['genelist', 'uniprotIds', 'GO_metadata', 
                               'GO_metadata_compartment', 'ECO_metadata_compartment']).dropna().keys())

        if kwList == set(['cytosol']):
            finalLocation = 'Cytosol'
        elif kwList == set(['cytosol', 'cytoplasm']):
            finalLocation = 'Cytosol'
        ############################################

        elif kwList == set(['cytosol', 'membrane']):
            finalLocation = 'Cytosol (IM-Associated)' 
        elif kwList == set(['cytosol', 'inner membrane', 'membrane']):
            finalLocation = 'Cytosol (IM-Associated)'
        elif kwList == set(['cytosol', 'cytoplasm', 'membrane']):
            finalLocation = 'Cytosol (IM-Associated)'
        elif kwList == set(['cytosol', 'cytoplasm', 'inner membrane', 'membrane']):
            finalLocation = 'Cytosol (IM-Associated)'
        elif kwList == set(['cytosol', 'cytoplasm', 'inner membrane', 'transmembrane', 'membrane']):
            finalLocation = 'Cytosol (IM-Associated)'
        ############################################
        elif kwList == set(['inner membrane', 'membrane']):
            finalLocation = 'Inner Membrane (No Area)'
        ############################################
        elif kwList == set(['periplasmic space', 'inner membrane', 'periplasm', 'membrane']):
            finalLocation = 'Periplasm (IM-Associated)'
        ############################################
        elif kwList == set(['periplasmic space']):
            finalLocation = 'Periplasm'
        elif kwList == set(['periplasmic space', 'periplasm']):
            finalLocation = 'Periplasm'
        ############################################
        elif kwList == set(['periplasmic space', 'outer membrane', 'membrane-bounded', 'periplasm', 'membrane']):
            finalLocation = 'Periplasm (OM-Associated)'
        ############################################
        elif kwList == set(['outer membrane', 'membrane']):
            finalLocation = 'Outer Membrane (No Area)'
        ############################################
        elif kwList == set(['extracellular']):
            finalLocation = 'Extracellular (secreted)'
        ############################################
        elif kwList ==set(['membrane']):
            finalLocation = 'Membrane Unknown (No Area)'
        ############################################
        elif 'pilus' in str(dforientation.loc[index].values):
            finalLocation = 'Periplasm'

        elif kwList ==set(['periplasmic space', 'periplasm', 'cytoplasm', 'cytosol']):
            finalLocation = 'Periplasm'
            
        elif kwList ==set(['periplasmic space', 'periplasm', 'cytosol']):
            finalLocation = 'Periplasm'

            
        

        else:
            print ('Could not determine compartment for ',c, index, kwList)
            
            indexmissing +=[index]
            c +=1
            continue

        dforientation.loc[index,'finalLocation'] = finalLocation
    return dforientation

def run_005F_addComparmentsManually(dforientation,
                                    manualOrientFile = op.join(qspaceDirs['Input_dir'],'005F-ManualCompartmentsForNonMembraneStructures.csv'),
                                    
                                   ):
    
    if op.exists(manualOrientFile):
        dfmanual_compartment  =pd.read_csv(manualOrientFile)
    else:
#         print "No Manual Compartment File Exists..."
        outfile = op.join(qspaceDirs['DataOutput_dir'],'005F-df_orientation.csv')
        dforientation.to_csv(outfile)
#         print "Saving....\n\t> {}".format(outfile)
        return dforientation
    
    for index, row in dfmanual_compartment.iterrows():
        if index in dforientation.index:
            dforientation.loc[index,'finalLocation'] = row.Manual
            
    outfile = op.join(qspaceDirs['DataOutput_dir'],'005F-df_orientation_after_manual_curation.csv')
    dforientation.to_csv(outfile)
    log.info("Saving Manual orientation of Proteins dataframe...\n\t{}".format(outfile))
    return dforientation


def run_005F_assignAALevelLocation(dforientation,
                                   dfalldata,
                                   data_int_orient,
                                   data_int_embed,
                                  ):
    
    
    compartment_dict = {} #non-membrane
    for structureId in tqdm(dforientation.index.tolist()):
        row = dforientation.loc[structureId]
        dfs = dfalldata[dfalldata.structureId == structureId]
        for index in dfs.index:
            compartment_dict.update({index : row.finalLocation})
    for index, value in compartment_dict.items():
        if 'associated' in str(value):
            value = value.replace('associated', 'Associated')
            compartment_dict.update({index  : value})
        if 'TransMembrane' == str(value):
            compartment_dict.update({index  : 'Unknown'})
        if str(value) =='nan':
            compartment_dict.update({index  : 'Unknown'})
            
    print ("Assigning AA Compartments.....\n----------------------")       
    aa_compartment_dict = {}
    for index, closest_leaf in data_int_orient.items(): #membrane
        if closest_leaf in [u'Cytoplasmic',u'Periplasmic_C', u'IM_1', u'IM_2',]:
            newLabel = 'IM-{}'.format(closest_leaf)
        elif closest_leaf in [u'Extracellular',u'Periplasmic_E', u'OM_1', u'OM_2',]:
            newLabel = 'OM-{}'.format(closest_leaf)
        else:
#             print (closest_leaf)
            continue
        if data_int_embed[index] == True:
            newLabel += '-Embedded'
        else:
            newLabel += '-Bulb'
        aa_compartment_dict.update({index : newLabel})
    aa_compartment_dict.update(compartment_dict) #both membrane and non-membrane updated here

    outfile = op.join(qspaceDirs['DataOutput_dir'], "005F-StructuralProteome_AA_compartment.json") 
    with open(outfile, 'w') as f:
        json.dump(aa_compartment_dict, f)
    log.info("Saving QSPACE Proteome - AA-level compartments...\n\t{}".format(outfile))
    
    
    
    ##################################################################################
    #get Protein-level orientation
    dfalldata['AA_compartment'] = pd.Series(aa_compartment_dict)
    protein_location = {}

    for structureId in tqdm(dfalldata.structureId.unique()):
        dfs =  dfalldata[dfalldata.structureId == structureId]
        dfs = dfs[dfs.structNum.isna() == False]
        protein_location.update({structureId:{'AA_compartment': dfs.AA_compartment.unique(), 'indexlist' : dfs.index.tolist()}})
    
    df = pd.DataFrame.from_dict(protein_location, orient = 'index')
    for index, row in df[df.AA_compartment.isna() == False].iterrows():
    #     print row.get('AA_compartment')

        if len(row.get('AA_compartment')) > 1:
            locs = []
            for aa_comp in row.get("AA_compartment"):
                if str(aa_comp)=='nan':
                    continue

                locs += [aa_comp[0:2]]
            locs = set(locs)

            if locs ==set(['OM', 'IM']):
                protein_loc = 'TransMembrane'
            elif locs == set(['IM']):
                protein_loc = 'Inner Membrane'
            elif locs == set(['OM']):
                protein_loc = 'Outer Membrane'

    #         print protein_loc
            df.loc[index, 'protein_compartment'] = protein_loc
        elif len(row.get('AA_compartment')) == 0 :
            df.loc[index, 'protein_compartment'] = 'Unknown'

        else:
    #         print row.get('AA_compartment')
            df.loc[index, 'protein_compartment'] = row.get('AA_compartment')[0]
    #         print row.get('AA_compartment')
    #     print row.location
    
    protein_final_compartment  = {}
    for index, row in tqdm(df.iterrows()):
        for glob_index in row.indexlist:
            protein_final_compartment.update({glob_index : row.protein_compartment})
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], "005F-StructuralProteome_Protein_compartment.json") 
    with open(outfile, 'w') as f:
        json.dump(protein_final_compartment, f)
    log.info("Saving QSPACE Proteome - Proteome-level compartments...\n\t{}".format(outfile))
    
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], "005F-Cell_compartments_protein_level.csv") 
    df.to_csv(outfile)
#     print "Saving....\n\t> {}".format(outfile)
    
    return protein_final_compartment, aa_compartment_dict,df
    
