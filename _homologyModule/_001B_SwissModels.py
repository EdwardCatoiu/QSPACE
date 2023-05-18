# -*- coding: utf-8 -*-

from .utils import *
log = logging.getLogger(__name__)
        


def run_001B(blattnerUniprotMapping,
          swissRepositoryFolder,
          force_realign= False,
         ):
    """
    ### Description
    The purpose of this notebook is to:
    - move the swiss files from the SWISS-MODEL/REPOSITORY into our appropriate folder
    - determine the real chains in each SWISS model
    - find missing alignment files (gene-chain-structure alignment)

    ### Input(s)
    1. SWISS Model Repository

    ### Output(s)

    1. Download SWISS MODEL Structure Files (.pdb)
    2. Determine real protein chains in SWISS models
    3. Expected needle alignment files 


    ### Folder(s)

    ```
    -------------------------InPuts------------------------------------------------
    SwissModelRepository <User Defined>
    └───INDEX.csv
    |       ...
    -------------------------OutPuts------------------------------------------------
    qspace
    └───data
    │   |   2. 001B-SWISS-model_chains_to_genes.json
    │   |   3. 001B-expected_alignment_files_SWISS.json
    |       ...
    //
    externalHardDrive <User Defined>
    └───structures
    |   |
    |   └───all_swiss
    |   |   |   1. {SWISS-MODEL-ID}.pdb
    |   |       ...
    ```
    """
    
    dfblattner_uniprot_mapping = pd.DataFrame()
    dfblattner_uniprot_mapping['uniprotId'] = pd.Series(blattnerUniprotMapping)
    dfblattner_uniprot_mapping.head()
    
    dfswiss = pd.read_csv(op.join(swissRepositoryFolder, 'INDEX.csv'), skiprows=6, sep ='\t')
    dfswiss.head()
    
    
    
    swiss_outfolder = qspaceDirs['swissStructuresDir']
    print ("Moving SWISS Models from the repository into...\n------------------------------------------------\n\t> {}".format(swiss_outfolder))
    for index, row in tqdm(dfswiss.iterrows()):
        uniprotId = row.get('UniProtKB_ac')
        if uniprotId not in blattnerUniprotMapping.values():
            continue

        [f1,f2,f3] = [uniprotId[0:2], uniprotId[2:4], uniprotId[4:]]

        swiss_infolder = op.join(swissRepositoryFolder,'{}/{}/{}/swissmodel/'.format(f1,f2,f3) )
        structure_id = '{}_{}_{}_{}.pdb'.format(uniprotId, row.get('from'), row.get('to'), row.get("template").split('.')[0] )
        outfilepath = op.join(swiss_outfolder,structure_id)
        if op.exists(outfilepath):
            continue

        filename = '{}_{}_{}_{}.pdb'.format(row.get('from'), row.get('to'), row.get("template"), row.get('coordinate_id'))
        infilepath = op.join(swiss_infolder, filename)


        if not op.exists(infilepath):
#             print filename, '\t', structure_id
            continue
        else:
            shutil.copy(src= infilepath,dst = outfilepath)
#     print "\t> Moved : {} SWISS models...".format(len(os.listdir(swiss_outfolder)))
    
    total_chains_swiss = get_struct_chains_swiss(structureFolder = qspaceDirs['swissStructuresDir'],
                                                 dfblattner_uniprot_mapping= dfblattner_uniprot_mapping)
    
    outfile = op.join(qspaceDirs['DataOutput_dir'],'001B-SWISS-model_chains_to_genes.json')
    log.info('Saving SWISS-chain-gene maps ... \n\t{}'.format(outfile))
    with open(outfile, 'w') as f :
        json.dump(total_chains_swiss,f)

    outfilepath = op.join(qspaceDirs['DataOutput_dir'],'001B-expected_alignment_files_SWISS.json' )
    find_missing_needle_files(total_chains_dict=total_chains_swiss, 
                              outfilepath = outfilepath,
                              blattnerUniprotMapping=blattnerUniprotMapping,
                              SequenceAlignmentDir= qspaceDirs['SequenceAlignmentDir'],
                              force_rerun= force_realign, 
                             )