# -*- coding: utf-8 -*-


import os
import os.path as op
from tqdm import tqdm_notebook as tqdm
import pandas as pd
import time
import numpy as np
import json
import shutil
import requests
import ast
from Bio import SeqIO

with open('qspace_directories.json','r') as f:
    qspaceDirs= json.load(f)
    
import logging
log = logging.getLogger(__name__)

##################################
import sys

# import qspace_directories as qspaceDirs
#import ssbio Functions
# sys.path.append(qspaceDirs.ssbioDir)
import ssbio
# sys.path = ['/usr/local/lib/python3.7/dist-packages/ssbio/'] + sys.path
# print (ssbio.__path__)
# print (sys.path)
from ssbio.protein.structure.structprop import StructProp
from ssbio.protein.sequence.properties.residues import one_to_three
from ssbio.protein.structure.utils.structureio import StructureIO
# from ssbio.structure.utils.structureio import StructureIO

from ssbio.protein.structure.utils import cleanpdb

#import qspace functions
# sys.path.append(qspaceDirs.FunctionsDir)
import download_alphafold_monomers as downloadAlphaFold
import check_structure
import alignment
import clean_itasser_string
import pseudo_structures as pseudo
#########################################################################

codon_table = {'UUU':'F','UUC':'F','UUA':'L','UUG':'L',
               'CUU':'L','CUC':'L','CUA':'L','CUG':'L',
               'AUU':'I','AUC':'I','AUA':'I','AUG':'M',
               'GUU':'V','GUC':'V','GUA':'V','GUG':'V',
               'UCU':'S','UCC':'S','UCA':'S','UCG':'S',
               'CCU':'P','CCC':'P','CCA':'P','CCG':'P',
               'ACU':'T','ACC':'T','ACA':'T','ACG':'T',
               'GCU':'A','GCC':'A','GCA':'A','GCG':'A',
               'UAU':'Y','UAC':'Y','UAA':'*','UAG':'*',
               'CAU':'H','CAC':'H','CAA':'Q','CAG':'Q',
               'AAU':'N','AAC':'N','AAA':'K','AAG':'K',
               'GAU':'D','GAC':'D','GAA':'E','GAG':'E',
               'UGU':'C','UGC':'C','UGA':'*','UGG':'W',
               'CGU':'R','CGC':'R','CGA':'R','CGG':'R',
               'AGU':'S','AGC':'S','AGA':'R','AGG':'R',
               'GGU':'G','GGC':'G','GGA':'G','GGG':'G',
               }
amino_acid_three = []
for aa in set(codon_table.values()):
    try:
        amino_acid_three += [one_to_three(aa)]
    except KeyError:
        continue
#############################################################################

alphafold_is_better = {'AAEA_ECOLI_model1_clean_residues_removed': [u'AF-P46482-F1-model_v2'],
                       'ADRA_ECOLI_model1_clean_residues_removed': [u'AF-P0AAP1-F1-model_v2'],
                       'ADRB_ECOLI_model1_clean_residues_removed': [u'AF-P76261-F1-model_v2'],
                       'ANSP_ECOLI_model1_clean_residues_removed': [u'AF-P77610-F1-model_v2'],
                       'ASMA_ECOLI_model1_clean_residues_removed': [u'AF-P28249-F1-model_v2'],
                       'CSGE_ECOLI_model1_clean_residues_removed': [u'AF-P0AE95-F1-model_v2'],
                       'CVRA_ECOLI_model1_clean_residues_removed': [u'AF-P76007-F1-model_v2'],
                       'CYSP_ECOLI_model1_clean_residues_removed': [u'AF-P16700-F1-model_v2'],
                       'DDPA_ECOLI_model1_clean_residues_removed': [u'AF-P76128-F1-model_v2'],
                       'FECR_ECOLI_model1_clean_residues_removed': [u'AF-P23485-F1-model_v2'],
                       'FLGA_ECOLI_model1_clean_residues_removed': [u'AF-P75933-F1-model_v2'],
                       'GFCD_ECOLI_model1_clean_residues_removed': [u'AF-P75882-F1-model_v2'],
                       'HEMX_ECOLI_model1_clean_residues_removed': [u'AF-P09127-F1-model_v2'],
                       'HOFC_ECOLI_model1_clean_residues_removed': [u'AF-P36646-F1-model_v2'],
                       'HYFH_ECOLI_model1_clean_residues_removed': [u'AF-P77423-F1-model_v2'],
                       'KEFB_ECOLI_model1_clean_residues_removed': [u'AF-P45522-F1-model_v2'],
                       'LLDP_ECOLI_model1_clean_residues_removed': [u'AF-P33231-F1-model_v2'],
                       'MBHM_ECOLI_model1_clean_residues_removed': [u'AF-P0ACE0-F1-model_v2'],
                       'MDTE_ECOLI_model1_clean_residues_removed': [u'AF-P37636-F1-model_v2'],
                       'MDTO_ECOLI_model1_clean_residues_removed': [u'AF-P32715-F1-model_v2'],
                       'MDTP_ECOLI_model1_clean_residues_removed': [u'AF-P32714-F1-model_v2'],
                       'MLTA_ECOLI_model1_clean_residues_removed': [u'AF-P0A935-F1-model_v2'],
                       'NARY_ECOLI_model1_clean_residues_removed': [u'AF-P19318-F1-model_v2'],
                       'OPGC_ECOLI_model1_clean_residues_removed': [u'AF-P75920-F1-model_v2'],
                       'PBP2_ECOLI_model1_clean_residues_removed': [u'AF-P0AD65-F1-model_v2'],
                       'PERM_ECOLI_model1_clean_residues_removed': [u'AF-P0AFI9-F1-model_v2'],
                       'PITB_ECOLI_model1_clean_residues_removed': [u'AF-P43676-F1-model_v2'],
                       'PSD_ECOLI_model1_clean_residues_removed': [u'AF-P0A8K1-F1-model_v2'],
                       'PTIBC_ECOLI_model1_clean_residues_removed': [u'AF-P24241-F1-model_v2'],
                       'RNFB_ECOLI_model1_clean_residues_removed': [u'AF-P77223-F1-model_v2'],
                       'RNFG_ECOLI_model1_clean_residues_removed': [u'AF-P77285-F1-model_v2'],
                       'SLP_ECOLI_model1_clean_residues_removed': [u'AF-P37194-F1-model_v2'],
                       'SLYB_ECOLI_model1_clean_residues_removed': [u'AF-P0A905-F1-model_v2'],
                       'SMP_ECOLI_model1_clean_residues_removed': [u'AF-P0AGC7-F1-model_v2'],
                       'UBIB_ECOLI_model1_clean_residues_removed': [u'AF-P0A6A0-F1-model_v2'],
                       'VIAA_ECOLI_model1_clean_residues_removed': [u'AF-P0ADN0-F1-model_v2'],
                       'WCAD_ECOLI_model1_clean_residues_removed': [u'AF-P71238-F1-model_v2'],
                       'WZYE_ECOLI_model1_clean_residues_removed': [u'AF-P27835-F1-model_v2'],
                       'YAHN_ECOLI_model1_clean_residues_removed': [u'AF-P75693-F1-model_v2'],
                       'YAIZ_ECOLI_model1_clean_residues_removed': [u'AF-P0AAQ0-F1-model_v2'],
                       'YAJI_ECOLI_model1_clean_residues_removed': [u'AF-P46122-F1-model_v2'],
                       'YBAY_ECOLI_model1_clean_residues_removed': [u'AF-P77717-F1-model_v2'],
                       'YBDJ_ECOLI_model1_clean_residues_removed': [u'AF-P77506-F1-model_v2'],
                       'YBFN_ECOLI_model1_clean_residues_removed': [u'AF-P75734-F1-model_v2'],
                       'YBHM_ECOLI_model1_clean_residues_removed': [u'AF-P75769-F1-model_v2'],
                       'YBHO_ECOLI_model1_clean_residues_removed': [u'AF-P0AA84-F1-model_v2'],
                       'YBIO_ECOLI_model1_clean_residues_removed': [u'AF-P75783-F1-model_v2'],
                       'YBJC_ECOLI_model1_clean_residues_removed': [u'AF-P46119-F1-model_v2'],
                       'YBJL_ECOLI_model1_clean_residues_removed': [u'AF-P60869-F1-model_v2'],
                       'YCDU_ECOLI_model1_clean_residues_removed': [u'AF-P75910-F1-model_v2'],
                       'YCEB_ECOLI_model1_clean_residues_removed': [u'AF-P0AB26-F1-model_v2'],
                       'YCGV_ECOLI_model1_clean_residues_removed': [u'AF-P76017-F1-model_v2'],
                       'YCHO_ECOLI_model1_clean_residues_removed': [u'AF-P39165-F1-model_v2'],
                       'YCIQ_ECOLI_model1_clean_residues_removed': [u'AF-P45848-F1-model_v2'],
                       'YDAM_ECOLI_model1_clean_residues_removed': [u'AF-P77302-F1-model_v2'],
                       'YDBH_ECOLI_model1_clean_residues_removed': [u'AF-P52645-F1-model_v2'],
                       'YDHJ_ECOLI_model1_clean_residues_removed': [u'AF-P76185-F1-model_v2'],
                       'YDHK_ECOLI_model1_clean_residues_removed': [u'AF-P76186-F1-model_v2'],
                       'YDJN_ECOLI_model1_clean_residues_removed': [u'AF-P77529-F1-model_v2'],
                       'YDJX_ECOLI_model1_clean_residues_removed': [u'AF-P76219-F1-model_v2'],
                       'YEAJ_ECOLI_model1_clean_residues_removed': [u'AF-P76237-F1-model_v2'],
                       'YEAY_ECOLI_model1_clean_residues_removed': [u'AF-P0AA91-F1-model_v2'],
                       'YFAZ_ECOLI_model1_clean_residues_removed': [u'AF-P76471-F1-model_v2'],
                       'YFBS_ECOLI_model1_clean_residues_removed': [u'AF-P0AFU2-F1-model_v2'],
                       'YFCC_ECOLI_model1_clean_residues_removed': [u'AF-P39263-F1-model_v2'],
                       'YGAW_ECOLI_model1_clean_residues_removed': [u'AF-P64550-F1-model_v2'],
                       'YGIM_ECOLI_model1_clean_residues_removed': [u'AF-P0ADT8-F1-model_v2'],
                       'YHDP_ECOLI_model1_clean_residues_removed': [u'AF-P46474-F1-model_v2'],
                       'YHFK_ECOLI_model1_clean_residues_removed': [u'AF-P45537-F1-model_v2'],
                       'YHFT_ECOLI_model1_clean_residues_removed': [u'AF-P45546-F1-model_v2'],
                       'YHGE_ECOLI_model1_clean_residues_removed': [u'AF-P45804-F1-model_v2'],
                       'YIAD_ECOLI_model1_clean_residues_removed': [u'AF-P37665-F1-model_v2'],
                       'YIAT_ECOLI_model1_clean_residues_removed': [u'AF-P37681-F1-model_v2'],
                       'YIDE_ECOLI_model1_clean_residues_removed': [u'AF-P60872-F1-model_v2'],
                       'YIDX_ECOLI_model1_clean_residues_removed': [u'AF-P0ADM6-F1-model_v2'],
                       'YIHO_ECOLI_model1_clean_residues_removed': [u'AF-P32136-F1-model_v2'],
                       'YJBH_ECOLI_model1_clean_residues_removed': [u'AF-P32689-F1-model_v2'],
                       'YMFE_ECOLI_model1_clean_residues_removed': [u'AF-P75968-F1-model_v2'],
                       'YNIB_ECOLI_model1_clean_residues_removed': [u'AF-P76208-F1-model_v2'],
                       'YOJI_ECOLI_model1_clean_residues_removed': [u'AF-P33941-F1-model_v2'],
                       'YPJA_ECOLI_model1_clean_residues_removed': [u'AF-P52143-F1-model_v2'],
                       'YQGA_ECOLI_model1_clean_residues_removed': [u'AF-Q46831-F1-model_v2'],
                       'YQHA_ECOLI_model1_clean_residues_removed': [u'AF-P67244-F1-model_v2'],
                       'YRAQ_ECOLI_model1_clean_residues_removed': [u'AF-P45468-F1-model_v2']
}

def get_struct_chains_alphafold(structureFolder,dfalphafold_monomer):
#     print '------------------------------------------------'
#     print "Ensuring all structure subunits are real chains..." 
    total_chains_dict = {}
    for structureFile in tqdm(os.listdir(structureFolder)):

        structureFilePath = op.join(structureFolder, structureFile)
#         print structureFile,structureFilePath
        dfu = dfalphafold_monomer[dfalphafold_monomer.alphafold_sfile == structureFilePath]
        if len(dfu) ==0:
            continue

        uniprotId = dfu.uniprotId.values[0]
        gene = dfu.first_valid_index()

        s = StructProp(ident=op.basename(structureFilePath).split('.pdb')[0] , structure_path=structureFilePath,file_type = 'pdb')
        try:
            p =s.parse_structure()
        except AttributeError:
            log.warn("Error parsing chains of structure file : {}".format(structureFile))
#             print structureFile
            continue
        real_chains = check_structure.find_real_chains(p)

        structure_id = op.basename(structureFilePath).split('.pdb')[0]
        if structure_id not in total_chains_dict:
            total_chains_dict.update({structure_id : {gene : real_chains}})
        elif gene not in total_chains_swiss[structure_id]:
            total_chains_dict[structure_id].update({gene : real_chains})
    return total_chains_dict

def get_structure_chains_itasser(structureFolder,itasser_metadata,blattner_uniprot_mapping,):
#     print '------------------------------------------------'
#     print 'Finding all real protein chains in ITASSER MODELS...'
    
    total_chains = {}
    nf = 0
    f= 0 

    for gene, uniprotId in tqdm(blattner_uniprot_mapping.items()):
        try:
            itasser_entry = itasser_metadata.loc[uniprotId,'Entry Name']
        except KeyError:
            nf +=1
#             print nf ,'\t', gene,'\t',uniprotId,'missing ITASSER model'
            log.warn("No ITASSER entry found : {}-{}".format(gene,uniprotId))

            continue

        structure_file = op.join(structureFolder, '{}_model1_clean_residues_removed.pdb'.format(itasser_entry))

        if not op.exists(structure_file):
            continue
        f +=1

        s = StructProp(ident=op.basename(structure_file).split('.pdb')[0] , structure_path=structure_file,file_type = 'pdb')
        p =s.parse_structure()
        real_chains = check_structure.find_real_chains(p)

        structure_id = op.basename(structure_file).split('.pdb')[0]
        if structure_id not in total_chains:
            total_chains.update({structure_id : {gene : real_chains}})
        elif gene not in total_chains[structure_id]:
            total_chains[structure_id].update({gene : real_chains})
    return total_chains
    
def get_struct_chains_swiss(structureFolder, dfblattner_uniprot_mapping):
    print ('Finding all real protein chains in SWISS MODELS...\n-----------------------------------')

    total_chains_swiss = {}
    for swiss_file in tqdm(os.listdir(structureFolder)):
        structure_id = swiss_file.split('.pdb')[0]
    #     if structure_id !='P0A9D4_2_262_1t3d':
    #         continue

        uniprotId = swiss_file.split('_')[0]  
        dfu = dfblattner_uniprot_mapping[dfblattner_uniprot_mapping.uniprotId == uniprotId]
        if len(dfu) == 0:
            continue

        if len(dfu) !=1:
#             print swiss_file
            continue
        gene = dfu.index.values[0]
        if len(dfu.index.unique()) > 1:
            raise ValueError

        structure_file = op.join(structureFolder, swiss_file)
        s = StructProp(ident=op.basename(swiss_file).split('.pdb')[0] , structure_path=structure_file,file_type = 'pdb')
#         p =s.parse_structure_old()
        p =s.parse_structure()
        
        real_chains = check_structure.find_real_chains(p)
        if structure_id not in total_chains_swiss:
            total_chains_swiss.update({structure_id : {gene : real_chains}})
        elif gene not in total_chains_swiss[structure_id]:
            total_chains_swiss[structure_id].update({gene : real_chains})
    #     if structure_id =='P0A9D4_2_262_1t3d':
    #         print structure_id, real_chains
    return total_chains_swiss

def find_missing_needle_files(total_chains_dict, 
                              outfilepath,
                              blattnerUniprotMapping,
                              SequenceAlignmentDir= qspaceDirs['SequenceAlignmentDir'],
                              force_rerun= False,
                              
                             ):
#     print '------------------------------------------------'
    print ('Finding missing needle alignments ...\n--------------------------------------')

    struct_error = []
    error_needle = []
    expected_files = {}
    for pdb_id in tqdm(total_chains_dict):
        missing_needle_files = alignment.find_expected_needle_files(pdb_id,
                                                                     total_chains=total_chains_dict,
                                                                     blattner_uniprot_mapping=blattnerUniprotMapping,
                                                                     structure_error=struct_error, 
                                                                     force_rerun = force_rerun,
                                                                     error_needle = error_needle,
                                                                     seq_alignment_folder=SequenceAlignmentDir,
                                                                    )
        if missing_needle_files !=[]:
            expected_files.update({pdb_id : missing_needle_files})


#     outfile = op.join(DataOutputDir,'001A-expected_alignment_files_ALPHAFOLD.json' )
#     print "{} Missing Needle alignments".format(len(expected_files))
#     print "\t> {}".format(outfilepath)

    with open(outfilepath,'w') as f:
        json.dump(expected_files, f)
    log.info('Saving Missing Needle Alignments...\n\t{}'.format(outfilepath))
