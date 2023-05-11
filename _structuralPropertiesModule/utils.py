import pandas as pd
from tqdm import tqdm_notebook as tqdm
import os
import os.path as op
import json
# import shutil
import sys
import numpy as np
import ast
import copy
# #import qspace directories
import sys
import logging

with open('qspace_directories.json','r') as f:
    qspaceDirs= json.load(f)
# #import qspace functions
sys.path.append(qspaceDirs['FunctionsDir'])
# import read_uniprot_text
# import uniprot_features
import protein_geometry
# # import prepare_structures_for_BFS as prepBFS
# # import oligomerization
# # import bfs_utils

import subprocess
import tempfile
from Bio import SeqIO


# import ssbio Functions
import ssbio
from ssbio.protein.structure.structprop import StructProp
from ssbio.protein.sequence.seqprop import SeqProp
from Bio.PDB.PDBExceptions import PDBConstructionException
from ssbio.protein.structure.utils.structureio import StructureIO


from Bio.PDB.Polypeptide import three_to_one, one_to_three
# import Bio.PDB.DSSP as DSSP
# import Bio.PDB.PDBParser as PDBParser
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

# from Bio.PDB.DSSP import dssp_dict_from_pdb_file
# from ssbio.protein.structure.properties.dssp import get_dssp_df_on_file
from ssbio.protein.structure.properties.dssp import secondary_structure_summary
from ssbio.protein.structure.properties.dssp import calc_surface_buried
from Bio.PDB.Polypeptide import aa1
from Bio.PDB.DSSP import residue_max_acc
from Bio.PDB.DSSP import dssp_dict_from_pdb_file


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
amino_acid_ids_three = []
for codon, aa_1 in codon_table.items():
    try:
        
        amino_acid_ids_three += [one_to_three(aa_1)]
    except KeyError:
        continue
amino_acid_ids_three = list(set(amino_acid_ids_three))

# import requests
# from bs4 import BeautifulSoup
# from Bio.PDB.Polypeptide import three_to_one, one_to_three
from Bio.PDB.PDBExceptions import PDBConstructionException

def get_fastaFileLocations(dfalldata,
                           alleleomeDir = qspaceDirs['AlleleomeSeqsDir'],
                           uniprotDir = qspaceDirs['UniprotSeqsDir'],
                          ):
    dfFastaLoc = pd.DataFrame()
    for gene in tqdm(dfalldata.gene.unique()):
        dfg = dfalldata[dfalldata.gene == gene]
        if len(dfg.geneSeqId.unique()) > 1:
#             print (g)
            raise KeyError

        geneSeqId = dfg.geneSeqId.unique().tolist()[0]
        if 'WT_consensus' in geneSeqId:
            fasta_file =  op.join(alleleomeDir, '{}_WT_consensus.fasta'.format(gene))
        else:
            fasta_file =  op.join(uniprotDir, '{}.fasta'.format(geneSeqId))

        dfFastaLoc.loc[gene,'fasta_id'] = geneSeqId
        dfFastaLoc.loc[gene,'fasta_location'] = fasta_file
    return dfFastaLoc
        
def find_sfile(structureId, stype):
    if stype =='PDB':
        sfile = op.join(qspaceDirs['bioassemblyStructuresDir'], '{}.cif'.format(structureId)) 
        if op.exists(sfile):
            return sfile
        
        sfile = op.join(qspaceDirs['cifStructuresDir'], '{}.cif'.format(structureId)) 
        if op.exists(sfile):
            return sfile
        
        sfile = op.join(qspaceDirs['pdbStructuresDir'], '{}.pdb'.format(structureId)) 
        if op.exists(sfile):
            return sfile
        return False
    if stype =='SWISS':
        sfile = op.join(qspaceDirs['swissStructuresDir'], '{}.pdb'.format(structureId)) 
        if op.exists(sfile):
            return sfile
        return False
    
    if stype =='ITASSER':
        sfile = op.join(qspaceDirs['itasserCleanStringRemovedStructuresDir'], '{}.pdb'.format(structureId)) 
        if op.exists(sfile):
            return sfile
        return False
    
    if stype =='ALPHAFOLD':
        sfile = op.join(qspaceDirs['alphafoldStructuresDir'], '{}.pdb'.format(structureId)) 
        if op.exists(sfile):
            return sfile
        return False
    if stype =='ALPHAFOLD_MULTIMER':
        alphaStructId = structureId.split('_Alpha')[0]
        
        resultFolder = op.join(qspaceDirs['alphafoldMultimerDir'],'{}.result/'.format(alphaStructId))
#         print resultFolder
        
        if not op.exists(resultFolder):
            return False
        
        for f in os.listdir(resultFolder):
            if alphaStructId in f and ".pdb" in f:
                if 'rank_1' in f or "rank_001" in f:
                    sfile = op.join(resultFolder, f)
                    return sfile
        return False
    
    
def pdbLine(new_chain, residue, atom, i):
    new_line = '{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}         {:>2s}{:2s}\n'.format('ATOM',
                                                                                                                  i,
                                                                                                                  atom.fullname,
                                                                                                                  atom.altloc,
                                                                                                                  residue.resname,
                                                                                                                  new_chain,
                                                                                                                  residue.id[1],
                                                                                                                  '',
                                                                                                                  float(atom.coord[0]),
                                                                                                                  float(atom.coord[1]),
                                                                                                                  float(atom.coord[2]),
                                                                                                                  atom.occupancy,
                                                                                                                  atom.bfactor,
                                                                                                                  atom.element,
                                                                                                                  ""
                                                                                                                 )
    return new_line
    


def write_pdb(s, p, outfile ,only_real_chains= True):

    with open(outfile ,'wb') as f:
        f.write('MODEL \n'.encode())

        i = 1

        used_chains = {}
        for new_chain in p.first_model.get_chains():
#             print '\t',new_chain
            if new_chain.id in ['_','.','-',' ']:
                continue

            for residue in new_chain.get_residues():
                if residue.resname not in amino_acid_ids_three and only_real_chains:
                    continue
    #                 print residue.id[1],
                for atom in residue.get_atoms():
        #             print chain.id, residue.id, atom.id
                    newline = pdbLine(new_chain.id, residue, atom, i)
                    f.write(newline.encode())
                    
                    i +=1

            f.write('TER \n'.encode())
        f.write('END \n'.encode())
    
    return outfile
#   
#   