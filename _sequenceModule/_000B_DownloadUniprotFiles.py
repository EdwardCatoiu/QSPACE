# -*- coding: utf-8 -*-
from .utils import *

log = logging.getLogger(__name__)

########################################################################################

def run_000B(blattnerUniprotMapping,
        force_rerun = False
       ):
    """
    ### Description
    The purpose of this notebook is to use the blattner-uniprot mappings generated from <u>Notebook #000A</u> to download the 
    UniProt.txt and UniProt.fasta files into ~/qspace/GEMPRO/sequences/uniprot/

    ### Input(s)
    1. {Blattner : UniProt ID} dictionary from #000A

    ### Output(s)

    1. {Uniprot ID}.txt file
    2. {Uniprot ID}.fasta file

    ### Folder(s)
    ```
    -------------------------Inputs-------------------------
    qspace 
    └───data
    │   |   000A-uniprot_blattner_mapping.json
    |   |   ...
    |
    -------------------------Outputs-------------------------
    └───GEMPRO
    │   └───sequences
    │       └───uniprot
    │           | UniprotID.txt
    │           | UniprotID.fasta
    │           | ...
    ```
    """
#     print ("<<<<<<<<<<<<<<<<<<<<<<<<#000B>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        
    print ('Downloading Uniprot .txt and .fasta files for {} genes\n----------------------------------------------------'.format(len(blattnerUniprotMapping)))
    log.info('Saving files to :\n\t{}'.format(qspaceDirs["UniprotSeqsDir"]))

#     print ('\t> {}'.format(qspaceDirs["UniprotSeqsDir"]))

    for blattnerId, uniprotId in tqdm(blattnerUniprotMapping.items()):
        download_uniprot_file(uniprot_id=uniprotId, filetype ='txt', outdir=qspaceDirs["UniprotSeqsDir"], force_rerun=force_rerun)
        download_uniprot_file(uniprot_id=uniprotId, filetype ='fasta', outdir=qspaceDirs["UniprotSeqsDir"], force_rerun=force_rerun)
