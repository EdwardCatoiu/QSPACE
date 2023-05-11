# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)


########################################################################################

def run_000C(alleleomeInputFolder,
          geneListInput = '000A-gene_list_4349.txt',
        force_rerun = False
       ):
    """
    ### Description
    The purpose of this notebook is to use write fasta files for genes in the WT alleleome.:

    E.A. Catoiu, P. Phaneuf, J.M. Monk, B.O. Palsson, Whole genome sequences from wild-type and laboratory evolved strains define the 
    alleleome and establish its hallmarks. PNAS. 120(15), 2023. https://doi.org/10.1073/pnas.221883512


    ### Input(s)
    1. Gene list of Blattner IDs >> 000A-gene_list_4349.txt
    2. Alleleome WT consensus sequences from the alleome project >> {geneId}_dfz.csv 

    ### Output(s)

    1. WT consensus AA sequences for each gene >> {geneId}_WT_consensus.fasta

    ### Folder(s)

    ```
    -------------------------Inputs------------------------------------------------
    Alleleome Project
    └───Alleleome_Data
    |   └───dfz
    │       |   {geneId}_dfz.csv
    |           ...
    //
    qspace 
    └───inputs
    │   |   000A-gene_list_4349.txt
    |   |   ...
    -------------------------Outputs------------------------------------------------
    |
    └───GEMPRO
    │   └───sequences
    │       └───alleleome
    │           | {geneId}_WT_consensus.fasta
    │           | ...
    ```
    """
#     print ("<<<<<<<<<<<<<<<<<<<<<<<<#000C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    print ("Downloading Alleleome WT sequences from list of geneIds...\n-------------------------------------------------------------")

    infile = op.join(qspaceDirs['Input_dir'], geneListInput)
    f = open(infile)
    lines = f.readlines()
#     print lines[0:5]
    genelist = []
    for gene in lines:
        genelist +=[gene.split('\n')[0].split(' #')[0]]
    f.close()
    log.info('Gene list input file :\n\t{}'.format(infile))
    
#     print ("> List of {} Genes: \t {}".format(len(genelist), infile))

    log.info('Alleleome Data location :\n\t{}'.format(alleleomeInputFolder))

#     print "\t> {}"
#     alleleomeInputFolder='/home/ecatoiu/Projects/Nature_Alleleome/Alleleome_Data/dfz'
    print ("\nMoving Files...")
#     print ("> Moving files to : \t{}".format(qspaceDirs.AlleleomeDir))
    for g in tqdm(genelist):

        outfile = op.join(qspaceDirs['AlleleomeSeqsDir'],"{}_WT_consensus.fasta".format(g) )    
        if op.exists(outfile) and not force_rerun:
            continue

        dfzFile = op.join(alleleomeInputFolder, "{}_dfz.csv".format(g))
        if not op.exists(dfzFile):
            log.warn('Missing alleleome file :\n{}'.format(dfzFile))

#             print ('missing file : {}'.format(dfzFile))
            continue

        dfz =  pd.read_csv(dfzFile, index_col=0)

        consensus_seq  = "".join(dfz[dfz.consensus_AAseq!= '-'].consensus_AAseq.tolist())
        consensus_seq = consensus_seq.split('*')[0]

        with open(outfile, "w") as f:
            f.write(">{}_consensus".format(g) + "\n" + consensus_seq + "\n")

    log.info("Alleleome Seqs Downloaded for {} Genes".format(len(os.listdir(qspaceDirs['AlleleomeSeqsDir']))))

