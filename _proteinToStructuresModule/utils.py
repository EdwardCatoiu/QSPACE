import os
import os.path as op
import json
from tqdm import tqdm_notebook as tqdm
import ast
import pandas as pd
import numpy as np
# import urllib
# import urllib2
# import requests
import copy

from Bio.PDB.Polypeptide import three_to_one,one_to_three
from Bio import SeqIO

#import qspace directories
import sys

with open('qspace_directories.json','r') as f:
    qspaceDirs= json.load(f)
    
import logging

#import qspace functions
# sys.path.append(qspaceDirs.FunctionsDir)
import prepare_structures_for_BFS as prepBFS
import oligomerization
import bfs_utils


#import ssbio Functions
# sys.path.append(qspaceDirs.ssbioDir)
# import ssbio
from ssbio.protein.sequence.seqprop import SeqProp
from ssbio.protein.structure.structprop import StructProp
import ssbio.protein.sequence.utils.alignment as ssbioaln


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
    
def map_repseq_resnums_to_struct_resnums(seq_resnums,repchain_resnums,struct_prop,chain_id):
    to_repchain_index = {}
    #             print structure_id
    for seq_res in seq_resnums:
#         print seq_res
        try:
            ix = repchain_resnums[int(seq_res) - 1] - 1
#             print ix
        except ValueError:
#             print ('{}, {}: INDEX WARNING: no equivalent residue found in structure sequence'.format(structure_id , seq_res))
            continue
        if np.isnan(ix):
            print ('{}, {}: no equivalent residue found in structure sequence'.format(structure_id , seq_res))
            pass
        else:
            to_repchain_index[seq_res] = int(ix)
    chain = struct_prop.chains.get_by_id(chain_id)
    repchain_structure_mapping = chain.seq_record.letter_annotations['structure_resnums']
#     print "repchain_structure_mapping"
#     print len(repchain_structure_mapping)
#     print repchain_structure_mapping
    to_structure_resnums = {}
#     print 'to_repchain_index'
#     print to_repchain_index
#     print '\n'
    for i, k in enumerate(seq_resnums):
#         print i,k
        if str(k) == 'nan':
            continue
        k  = int(k)
#                     aa1 = seq_residue_AA1[i]
#                     try:
#                         aa3 = one_to_three(aa1)
#                     except KeyError:
#                         aa3 = 'XXX'
        v = np.nan
        rn = ('' , np.nan, '')
        try:
#             print k,
            v = to_repchain_index[k]
#             print v,
#             print v,
            rn = repchain_structure_mapping[v]
#             print rn
        except KeyError:
            pass
        if rn[1] == float('Inf'):
            #print ('{}, {}: structure file does not contain coordinates for this residue'.format(a_structure.id, k))
            pass
        else:
            to_structure_resnums[k] = rn
            
    return to_structure_resnums
            
        