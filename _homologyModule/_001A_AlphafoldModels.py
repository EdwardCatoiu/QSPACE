# -*- coding: utf-8 -*-
from .utils import *
# import utils

log = logging.getLogger(__name__)

        
def run_001A(blattnerUniprotMapping,
          force_redownload= False,
          force_realign= False,
          alphafoldStructuresDir= False,
          alphafoldMetricsDir= False
         ):
    """
    ### Description
    The purpose of this notebook is to:
    - download the AlphaFold (https://alphafold.ebi.ac.uk/) monomeric structural files and model metrics using the UniprotIDs 
    - determine the quality of each alphafold model (b-factor avg)
    - find missing alignment files (gene-chain-structure alignment)

    ### Input(s)
    1. {Blattner : UniProt ID} dictionary from #000A

    ### Output(s)
    1. Download Alphafold Structure Files (.pdb)
    2. Download Alphafold Quality Metrics (dictionary)
    3. Quality of alphafold Monomers (dictionary)
    4. Quality of alphafold Monomers (.csv)
    5. Figure of alphafold Quality 
    6. Expected needle alignment files 

    ### Folder(s)
    ```
    -------------------------InPuts------------------------------------------------
    qspace 
    └───data
    │   |   1. 000A-uniprot_blattner_mapping.json
    |       ...
    -------------------------OutPuts------------------------------------------------
    qspace
    └───data
    │   |   3. 001A-quality_of_alphafold_monomers.json
    │   |   4. 001A-quality_of_alphafold_monomers.csv
    │   |   6. 001A-expected_alignment_files_ALPHAFOLD.json
    |       ...
    └───figures
    │   |   5. 001A-alphafold_structure_quality.png
            ...
    //
    externalHardDrive <User Defined>
    └───structures
    |   |
    |   └───all_alphafold
    |   |   |   1. AF-{uniprotId}-F1-model_v2.pdb
    |   |       ...
    |   └───all_alphafold-metrics
    |       |   2. AF-{uniprotId}-F1-predicted_aligned_error_v2.json
    |           ...
    ```
    """
    
    if not alphafoldStructuresDir:
        alphafoldStructuresDir = qspaceDirs['alphafoldStructuresDir']
#         print 'Using Default Directories....'
#         print '------------------------------------------------'
#     print 'Downloading Alphafold models to...'
#     print '\t> {}'.format(alphafoldStructuresDir)
    if not alphafoldMetricsDir:
        alphafoldMetricsDir = qspaceDirs['alphafoldMetricsDir']
#         print '------------------------------------------------'
#     print 'Downloading Alphafold metrics to...'
#     print '\t> {}'.format(alphafoldMetricsDir)

#     print alphafoldStructuresDir
#     print alphafoldMetricsDir

    for folder in [alphafoldStructuresDir, alphafoldMetricsDir]:
        if not op.exists(folder):
            os.mkdir(folder)

    print ("Downloading Alphafold Models and Metrics...\n----------------------------------------------------------")
    for blattnerId, uniprotId in tqdm(blattnerUniprotMapping.items()):
        #download the PDB file

        structure_url = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v2.pdb".format(uniprotId)
        outfolder=alphafoldStructuresDir
        if not op.exists(op.join(outfolder, op.basename(structure_url))):
            downloadAlphaFold.download_alphafold_structure(url= structure_url,outfolder=outfolder,force_rerun = force_redownload)
            time.sleep(0.5)

        #download the quality metrics
        metrics_url   = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-predicted_aligned_error_v2.json".format(uniprotId)
        outfolder=alphafoldMetricsDir
        if not op.exists(op.join(outfolder, op.basename(metrics_url))):
            downloadAlphaFold.download_alphafold_metrics(url= metrics_url,outfolder=outfolder,force_rerun = force_redownload)
            time.sleep(0.5)
            
#     print '------------------------------------------------'

    print ('Checking Quality of Alphafold Structure Database Models...\n--------------------------------------------')
    dfalphafold_monomer = pd.DataFrame()
    bfactor_dict = {}
    uniprot_dict = {}
    sfilePath_dict = {}


    for blattnerId, uniprotId in tqdm(blattnerUniprotMapping.items()):
        structure_url = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v2.pdb".format(uniprotId)
        sfile_downloaded = op.basename(structure_url)
        sfile_path = op.join(alphafoldStructuresDir ,sfile_downloaded)
        if not op.exists(sfile_path):
            continue
    #     dfalphafold_monomer.loc[blattnerId,'uniprotId'] = uniprotId
    #     dfalphafold_monomer.loc[blattnerId,'alphafold_sfile'] = sfile_path

        uniprot_dict.update({blattnerId : uniprotId})
        sfilePath_dict.update({blattnerId : sfile_path})
        try:
            s = StructureIO(structure_file = sfile_path)
            s.first_model['A'].get_residues()
        except (IOError, AttributeError) as e:
#             print index, uniprotId, sfile
            continue

        #get the prediction confidence per reside (CA-atom bfactors)
        bfactor_list = []
        for residue in s.first_model['A'].get_residues():
            if residue.resname not in amino_acid_three:
                continue
            bfactor = residue['CA'].bfactor
            bfactor_list += [bfactor]
        bfactor_dict.update({blattnerId : np.mean(bfactor_list)}) 

    dfalphafold_monomer['uniprotId'] = pd.Series(uniprot_dict)
    dfalphafold_monomer['bfactor_avg'] = pd.Series(bfactor_dict)
    dfalphafold_monomer['alphafold_sfile'] = pd.Series(sfilePath_dict)
    
    
    
    DataOutputDir = qspaceDirs['DataOutput_dir']

    quality_of_structure = {}
    for index, row in dfalphafold_monomer[dfalphafold_monomer.bfactor_avg.isna()== False].iterrows():
        s_id = op.basename(dfalphafold_monomer.loc[index,'alphafold_sfile']).split('.pdb')[0]
        quality_of_structure.update({s_id : {index :{'A' : [ row.bfactor_avg > 50. , row.bfactor_avg]}}})

    outfile = op.join(DataOutputDir, '001A-quality_of_alphafold_monomers.json')
    log.info('Saving Quality of Alphafold Monomers... \n\t{}'.format(outfile))

    with open(outfile, 'w') as f:
        json.dump(quality_of_structure, f)

    outfile = op.join(DataOutputDir, '001A-quality_of_alphafold_monomers.csv')
#     print '\t> {}'.format(outfile)
    dfalphafold_monomer = dfalphafold_monomer[dfalphafold_monomer.bfactor_avg.isna()== False]
    dfalphafold_monomer.to_csv(outfile)
    log.info('Saving Quality of Alphafold Monomers... \n\t{}'.format(outfile))

    
    total_chains_dict= get_struct_chains_alphafold(structureFolder =alphafoldStructuresDir,
                                                   dfalphafold_monomer=dfalphafold_monomer)
    
    outfilepath = op.join(DataOutputDir,'001A-expected_alignment_files_ALPHAFOLD.json' )
    find_missing_needle_files(total_chains_dict=total_chains_dict, 
                              outfilepath = outfilepath,
                              blattnerUniprotMapping=blattnerUniprotMapping,
                              SequenceAlignmentDir= qspaceDirs['SequenceAlignmentDir'],
                              force_rerun= force_realign, 
                             )
    
    return dfalphafold_monomer