# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)



def parse_MSMS(msms_file, df):
    dfMSMS = pd.read_csv(msms_file)
    res_depth  = dfMSMS['res_depth']
    ca_depth  = dfMSMS['ca_depth']
    msms_info = map(lambda index, res_depth_info, ca_depth_info : ({index :res_depth_info}, {index:ca_depth_info}), df.index.values, res_depth.values, ca_depth.values)
    return msms_info
    

def run_006C_MSMS(dfalldata,
                  MSMSoutfolder = qspaceDirs['msmsResultsDir'],
                  outPDBfolder= qspaceDirs['cifToPdbFolder'],
                  force_rerun = False,
                  force_rewrite_pdb = False
                 ):
    
    
    print ("Getting structure file locations...")
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

    
    for structure_id in tqdm(dfalldata.structureId.unique()):
    #     if structure_id in [ 'P0AER0_6_259_1ldf','P30125_1_363_1cm7','1zn1-assembly1']:
    #         continue
        outfilename = '{}_msms.df'.format(structure_id)
        outfile = op.join(MSMSoutfolder, outfilename)
        if op.exists(outfile) and not force_rerun:
            continue

        file_location = dfstructure_file_locations.loc[structure_id, "sfile"]

        outfilename =  op.basename(file_location).split('.')[0] + '_msms.df'
        outfile = op.join(MSMSoutfolder, outfilename)
        if op.exists(outfile) and not force_rerun:
            continue

        print (structure_id)

        s = StructProp(ident=structure_id, structure_file=file_location)
        parsed = s.parse_structure()
    #     break
        if '.cif' in file_location:

            try:
                outfilePDB = op.join(outPDBfolder, structure_id + '.pdb') 
                if not op.exists(outfilePDB) or force_rewrite_pdb:
                    outfilePDB = write_pdb(s = s, p= parsed, outfile = outfilePDB)

                s = StructProp(ident=structure_id, structure_file=outfilePDB)
                parsed = s.parse_structure()
            except (PDBConstructionException,ValueError) as e:
                print ('MSMS error : {}'.format(structure_id))
                continue


        print (s.structure_path)
        msms_results = ssbio.structure.properties.msms.get_msms_df(model=parsed.first_model,
                                                                       pdb_file=s.structure_path,
                                                                       outdir=MSMSoutfolder,
                                                                       force_rerun=force_rerun)
        
    return dfstructure_file_locations
    
def run_006C_AAlevelMSMSannotation(dfalldata,dfstructure_file_locations):
    res_depth = {}
    ca_depth = {}
    msms_errorlist = []

    cplx_list = dfalldata.cplx.unique()
    for cplxId in tqdm(cplx_list):
        dfc = dfalldata[dfalldata.cplx ==cplxId]

        for structureId in dfc.structureId.unique():
            dfs = dfc[dfc.structureId ==structureId]
            dfs = dfs[dfs.structNum.isna() == False]
            dfs= dfs[dfs.structNum != np.inf]


            file_location = dfstructure_file_locations.loc[structureId, "sfile"]

            basename  = op.basename(file_location).split('.')[0]
            msms_file = op.join(MSMSoutfolder , '{}_msms.df'.format(basename))

            if not op.exists(msms_file):
                print ('error no msms file >>>' , structureId)
                msms_errorlist.append(structureId)
                continue

            msms_info = parse_MSMS(msms_file, df=dfs)
            for (r, c) in msms_info:
                res_depth.update(r)
                ca_depth.update(c)
                
    outfile = op.join(qspaceDirs['DataOutput_dir'],'006C-StructuralProteome_MSMS_ca_depth.json')    
    with open(outfile, 'w') as f:
        json.dump(ca_depth,f )
    log.info("Saving QSPACE Proteome - MS/MS CaDepth...\n\t{}".format(outfile))
    
    
    outfile = op.join(qspaceDirs['DataOutput_dir'],'006C-StructuralProteome_MSMS_res_depth.json')    
    with open(outfile, 'w') as f:
        json.dump(res_depth,f )
    log.info("Saving QSPACE Proteome - MS/MS ResDepth...\n\t{}".format(outfile))
    
    return ca_depth,res_depth 