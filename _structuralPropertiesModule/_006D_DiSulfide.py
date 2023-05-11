# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)



def read_SS_bridges(disulfide_bridges, s_id):
    i = 0
    disulfide_dict = {}
    for chain in disulfide_bridges:
        for (cys1,cys2) in disulfide_bridges[chain]:
            resnum1 = cys1[1]
            resnum2 = cys2[1]

            chain_res1 = "{}_{}".format(chain, resnum1)
            chain_res2 = "{}_{}".format(chain, resnum2)

            bond_id = "{}&{}".format(chain_res1,chain_res2)
            disulfide_dict.update({i: bond_id})
            i +=1
    disulfide_dict = {s_id : disulfide_dict}    
    return disulfide_dict

def run_006D_ssbioDisulfide(dfalldata,disulfideResultsFolder = qspaceDirs['disulfideDir'],threshold= 3.0):
    slist = dfalldata.structureId.unique().tolist()
    # slist.reverse()
    print ("Getting structure file locations...\n------------------------------------")
    structure_dict = {}

    for structure in dfalldata.structureId.unique():
        if 'assembly' in structure:
            stype = 'PDB'
        elif 'AlphaMulti' in structure:
            stype = 'ALPHAFOLD_MULTIMER'
        elif 'AF-' in structure and 'model' in structure:
            stype= 'ALPHAFOLD'
        elif 'ECOLI' in structure and 'clean_residues' in structure:
            stype= 'ITASSER'
        else:
            stype= 'SWISS'

        structure_dict.update({structure : stype})

    dfstructure_file_locations = pd.DataFrame()
    dfstructure_file_locations['stype'] = pd.Series(structure_dict)
    for index, row in dfstructure_file_locations.iterrows():
        sfile = find_sfile(index, row.stype)
        dfstructure_file_locations.loc[index, 'sfile'] =  sfile
    dfstructure_file_locations.head()
    
    print ("Finding Disulfide bridges ...\n------------------------------------")
    for structure_id in tqdm(slist):
    #     print structure_id
        file_location = dfstructure_file_locations.loc[structure_id,"sfile"]
        
        
        s = StructProp(ident=structure_id, structure_path=file_location, file_type = op.basename(file_location).split('.')[-1])
        parsed = s.parse_structure()   
        model = parsed.first_model
        disulfide_bridges = ssbio.protein.structure.properties.residues.search_ss_bonds(model,threshold=threshold)    
        if disulfide_bridges =={}:
            with open(op.join(disulfideResultsFolder, "{}_disulfide_bridges.json".format(structure_id)), 'w') as f:
                json.dump(disulfide_bridges,f)      
            continue

        else:
            disulfied_info = read_SS_bridges(disulfide_bridges, structure_id)
            with open(op.join(disulfideResultsFolder, "{}_disulfide_bridges.json".format(structure_id)), 'w') as f:
                json.dump(disulfied_info,f) 
            log.info("{} S-S bridges in {}".format( len(disulfied_info[structure_id]),structure_id))
            continue
    
    print ('AA-level annotation of Disulfide bridges ...\n------------------------------------')
    # disulfide_index = []
    global_disulfide_dict = {}
    for fileName in tqdm(os.listdir(disulfideResultsFolder)):
        filePath = op.join(disulfideResultsFolder, fileName)

        with open(filePath, 'r') as f:
            disulfide_dict = json.load(f)

            for structureId , DSinfo in disulfide_dict.items():
                if structureId not in dfalldata.structureId.unique():
                    continue
                for DSnum, DSbridge_string in DSinfo.items():
                    [cr1, cr2] = DSbridge_string.split('&')


                    dfs = dfalldata[dfalldata.structureId ==structureId]
                    if len(dfs) == 0:
#                         print structureId
                        continue
                    df_disulfide = dfs[dfs.sfileChain_Residue.isin([cr1, cr2])].index.tolist()
                    for index in df_disulfide:
                        global_disulfide_dict.update({int(index) : 1})
                 
    outfile = op.join(qspaceDirs['DataOutput_dir'],'006D-StructuralProteome-disulfide_bridges.json')    
    with open(outfile, 'w') as f:
        json.dump(global_disulfide_dict,f )
    log.info("Saving QSPACE Proteome - S-S bridges...\n\t{}".format(outfile))
    return global_disulfide_dict