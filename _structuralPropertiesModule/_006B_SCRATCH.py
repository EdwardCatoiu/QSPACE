# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)


def parse_SCRATCH(scratch_file, 
                  acc_type,
                  df,
                  scratchOutfolder = qspaceDirs['scratchResultsDir'],
                  ):
    

    scratchFilePath = op.join(scratchOutfolder,scratch_file)
    df_SCRATCH = pd.read_csv(scratchFilePath)

    scratch = df_SCRATCH.loc[0].values[0]
    if acc_type == 'acc20':
        scratch = scratch.split(' ')
        scratch = map(int, scratch)
    else:
        scratch = list(scratch)
        
#     print acc_type, len(scratch)
    scratch_dict = {}
    scratch_info = map(lambda index, info : {int(index) : info}, df.index.values, scratch)
    for index in scratch_info:
        scratch_dict.update(index)
    return scratch_dict
    
    
def run_006B_SCRATCH(dfalldata,
                     path_to_scratch,
                     num_cores = 4,
                     dfFastaLoc = pd.DataFrame(),
                     scratchFolderNEW = qspaceDirs['scratchResultsDir'],
                     force_rerun= False):
    
    if len(dfFastaLoc) ==0:
        dfFastaLoc = get_fastaFileLocations(dfalldata,
                                          alleleomeDir = qspaceDirs['AlleleomeSeqsDir'],
                                          uniprotDir = qspaceDirs['UniprotSeqsDir'],
                                         )
    
    
    print ("Running SCRATCH.....")
    indexlist = dfFastaLoc.index.tolist()
    for index in tqdm(indexlist):
        found_file = True
        row = dfFastaLoc.loc[index]


        if 'consensus' in row.fasta_id:
            fasta_id = '{}_consensus_allele'.format(index)
        else:
            fasta_id = '{}_{}'.format(index, row.fasta_id)

        fasta_file = dfFastaLoc.loc[index, 'fasta_location']
        outname = op.join(scratchFolderNEW, fasta_id )

        for ext in ['.ss','.ss8','.acc','.acc20']:
            if not op.exists(outname+ ext):
                print (outname)
                found_file = False
                break

        if not found_file or force_rerun:
            try:
                print (fasta_id)
                shell_command='{} {} {} {}'.format(path_to_scratch, fasta_file, outname, num_cores)
                ssbio.utils.command_runner(shell_command, force_rerun_flag=force_rerun, outfile_checker= outname)
            except ssbio.utils.subprocess.CalledProcessError:
                log.warn("Unable to run SCRATCH for {}".format(fasta_id))

                print ('ERROR^')

                
def run_006B_AAlevelSCRATCHannotation(dfalldata,
                                      scratchOutfolder = qspaceDirs['scratchResultsDir'],
                                     ):
    print ('AA-level SCRATCH annotation...')

    scratch_global_info = {}
    for cplxId in tqdm(dfalldata.cplx.unique()):
        dfc = dfalldata[dfalldata.cplx == cplxId]
        for structureId in dfc.structureId.unique():

            dfs = dfc[dfc.structureId == structureId]
            for geneId in dfs.gene.unique():
                dfg = dfs[dfs.gene == geneId]

                for chainId in dfg.structChainId_mod.unique():
                    dfgc = dfg[dfg.structChainId_mod == chainId]
                    uniprotId = dfgc.geneSeqId.values[0]

                    if 'consensus' in uniprotId:
                        scratchBaseName = '{}_consensus_allele'.format(geneId)
                    else:
                        scratchBaseName = '{}_{}'.format(geneId, uniprotId)



    #                 print geneId, uniprotId, chainId, scratchBaseName

                    for acc_type in ['ss','ss8','acc','acc20']:
                        scratchFileName = "{}.{}".format(scratchBaseName, acc_type)
                        scratch_info = parse_SCRATCH(scratch_file=scratchFileName,
                                                     acc_type=acc_type, 
                                                     df = dfgc,
                                                     scratchOutfolder = scratchOutfolder)
                        if acc_type in scratch_global_info:
                            scratch_global_info[acc_type].update(scratch_info)
                        else:
                            scratch_global_info.update({acc_type : scratch_info})
    
    outfile = op.join(qspaceDirs['DataOutput_dir'],'006B-StructuralProteome_SCRATCH.json')    
    with open(outfile, 'w') as f:
        json.dump(scratch_global_info,f )
#     log.info("Saving QSPACE AA-level Scratch results...\n\t{}".format(outfile))
    log.info("Saving QSPACE Proteome - AA-level Scratch results...\n\t{}".format(outfile))

    return scratch_global_info
