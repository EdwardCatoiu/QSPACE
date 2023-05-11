import pandas as pd
from tqdm import tqdm_notebook as tqdm
import os
import os.path as op
import json
import shutil
import sys
import numpy as np
#import qspace directories
import ast
import sys
# sys.path.append('../')
# import qspace_directories as qspaceDirs
import copy
# import qspace functions

# import prepare_structures_for_BFS as prepBFS
# import oligomerization
# import bfs_utils

with open('qspace_directories.json','r') as f:
    qspaceDirs= json.load(f)
    
sys.path.append(qspaceDirs['FunctionsDir'])
import read_uniprot_text
import uniprot_features
import protein_geometry

import logging
log = logging.getLogger(__name__)

#import ssbio Functions
# sys.path.append(qspaceDirs.ssbioDir)
import ssbio
from ssbio.protein.structure.structprop import StructProp
from ssbio.protein.sequence.seqprop import SeqProp

###### webserver OPM
from selenium.webdriver.chrome.options import Options
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
import time

import requests
from bs4 import BeautifulSoup
from Bio.PDB.Polypeptide import three_to_one, one_to_three
from Bio import SeqIO

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


def get_membrane_from_GO( dfrepseq, dfall,sequenceFolder=qspaceDirs['UniprotSeqsDir']):
    
    membraneGenesGO = []

    GO_terms_data = {}

    for index , row in tqdm(dfrepseq.iterrows()):
        txt_file = op.join(sequenceFolder, '{}.txt'.format(row.UniProtId))
        if not op.exists(txt_file):
            print (index, row.UniProtId,'uniprot text file DNE')
            continue
        p, f, c= read_uniprot_text.find_GO(txt_file)

        data = str(p).lower() + str(f).lower() +str(c).lower()

        GO_terms_data.update({index : data})

        ismembrane = False
        for keyword in ['membrane','transport','abc','wall','periplasm']:
            if keyword in data:
                membraneGenesGO +=[index]
                ismembrane = True
                break
    #     if not ismembrane:
    #         print p,f,c
    membraneStructuresGO = dfall[dfall.gene.isin(membraneGenesGO)].structureId.unique()

    outfile = op.join(qspaceDirs['DataOutput_dir'], "005A-membrane_kw_GO_terms.json")
    with open(outfile, 'w') as f:
        json.dump(GO_terms_data, f)
#     print "Saving....\n\t> {}".format(outfile)
    return set(membraneGenesGO), set(membraneStructuresGO),GO_terms_data 


def get_membrane_from_iML1515( dfiml1515, dfall):
    
    iml1515_data = {}
    for index, gene in tqdm(dfiml1515.m_gene.items()):
        iml1515_data.update({gene : dfiml1515.loc[index,'m_subsystem']})

    membraneGenesIml = []
    for index, row in dfiml1515.iterrows():
        for keyword in ['membrane','transport','abc','wall','periplasm','envelope','murein']:
            if keyword in row.m_subsystem.lower():
                #only add the gene if it is part of our query (dfrepseq)
#                 if row.m_gene in dfrepseq.index:
                membraneGenesIml += [row.m_gene]
                break

    membraneStructuresIml = dfall[dfall.gene.isin(membraneGenesIml)].structureId.unique()
#     print len(set(membraneGenesIml)), len(set(membraneStructuresIml))
    

    outfile = op.join(qspaceDirs['DataOutput_dir'], "005A-membrane_kw_IML1515_terms.json")
    with open(outfile, 'w') as f:
        json.dump(iml1515_data, f)
    log.info("Saving membrane genes identified from iML1515 keyword search...\n\t{}".format(outfile))
    return set(membraneGenesIml), set(membraneStructuresIml) 




def get_membrane_from_UniProt(dfrepseq,dfall,sequenceFolder=qspaceDirs['UniprotSeqsDir'],checked=[]):
    
    uniprot_mem_data = {}
    
    membraneGenesUNI = []
    for index in tqdm(dfrepseq.index.tolist()):
        if index in checked:
            continue
        row = dfrepseq.loc[index] 

        txt_file = op.join(sequenceFolder, '{}.txt'.format(row.UniProtId))
        faa_file = op.join(sequenceFolder, '{}.fasta'.format(row.UniProtId))

        if not op.exists(txt_file):
            print (index, row.UniProtId, 'uniprot text file DNE')
            continue
        if not op.exists(faa_file):
            print (index, row.UniProtId , 'uniprot fasta file DNE')
            continue

        feature_df = feature_df = uniprot_features.make_nextGenome_feature_dataframe(index,
                                                                                     uniprot_text_path=txt_file, 
                                                                                     uniprot_fasta_path=faa_file)
        isMembrane = False
        for feature in ['TOPO_DOM','INTRAMEM','TRANSMEM']:
            if feature in feature_df.index.tolist():
                membraneGenesUNI +=[index]
                isMembrane = True
                uniprot_mem_data.update({index : feature})
                break
    #     if not isMembrane:
    #         print feature_df.index.tolist()
        checked += [index]
    #     break
    membraneStructuresUNI = dfall[dfall.gene.isin(membraneGenesUNI)].structureId.unique()
    
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], "005A-membrane_kw_UniProt_terms.json")
    with open(outfile, 'w') as f:
        json.dump(uniprot_mem_data, f)
    log.info("Saving membrane structures identified from uniprot keyword search...\n\t{}".format(outfile))

    return set(membraneGenesUNI), set(membraneStructuresUNI)

def get_membrane_from_Ecocyc(dfall,
                            dfecocyc_loc = op.join(qspaceDirs['Input_dir'], "005A-EcocycSmartTable-All-genes-of-E.-coli-K-12-substr.-MG1655.txt"),                            
                            ):
    
    if type(dfecocyc_loc) == str:
        dfecocyc_loc = pd.read_csv(dfecocyc_loc, sep = '\t')
    dfecocyc_loc['gene'] = dfecocyc_loc['Accession-1']
    dfecocyc_loc = dfecocyc_loc[['Gene Name','gene','Locations']]
    dfecocyc_loc = dfecocyc_loc[dfecocyc_loc.gene.isna()==False]
    
    ECOCYC_loc = {}
    for index, row in dfecocyc_loc.iterrows():
        ECOCYC_loc.update({row.gene : row.Locations})
        
    potential_ecocyc = {}
    for index, row in tqdm(dfecocyc_loc.iterrows()):
        for kw in ['membrane','periplas','transport','secretion','extracellular','wall']:
            if kw in str(row.values).lower():
                potential_ecocyc.update({index : row.Locations})
                break
    dfecocyc_loc['potentialMembrane'] = pd.Series(potential_ecocyc)
    
    
    genes_ecocyc = list(dfecocyc_loc[dfecocyc_loc.potentialMembrane.isna() == False].gene.unique())
    structures_ecocyc = list(dfall[dfall.gene.isin(genes_ecocyc)].structureId.unique())
    
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], "005A-membrane_kw_Ecocyc_terms.json")
    with open(outfile, 'w') as f:
        json.dump(potential_ecocyc, f)
    log.info("Saving membrane structures identified from Ecocyc keyword search...\n\t{}".format(outfile))
    
    return set(genes_ecocyc), set(structures_ecocyc)


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
    
    
    
def get_pdb_coordinates(s_id, sfile):
    dfprotein_coords= pd.DataFrame()
    s = StructProp(ident = s_id, structure_path=sfile,file_type = op.basename(sfile).split('.')[-1])
    p = s.parse_structure()
    for c in p.first_model.get_chains():
        for res in c.get_residues():
            resname  = res.resname
            if resname in amino_acid_ids_three:
    #             print res.id, res['CA'].coord
                try:
                    dfprotein_coords.loc["{}_{}".format(c.id, res.id[1]), 'x'] = res['CA'].coord[0]
                    dfprotein_coords.loc["{}_{}".format(c.id, res.id[1]), 'y'] = res['CA'].coord[1]
                    dfprotein_coords.loc["{}_{}".format(c.id, res.id[1]), 'z'] = res['CA'].coord[2]
                except KeyError:
                    log.warn("Error residue {}-{} of {} has no coordinates......".format(c.id, res.id[1], s_id ))

                    continue
    return dfprotein_coords

def get_coordinates_of_opm_membrane(sfile, num_residues_to_consider  = 50):
    s = StructProp(ident = op.basename(sfile).split('.pdb')[0], structure_path=sfile, file_type = op.basename(sfile).split('.')[-1])
    p = s.parse_structure()
    N_coords = {}
    O_coords = {}
    resnum = 0 
    for c in p.first_model.get_chains():
        for res in c.get_residues():
            resname  = res.resname

            if resname in amino_acid_ids_three:
                continue
            if resname != 'DUM':
                continue
            if resnum > num_residues_to_consider:
                break
            
            
            try:
                N_coords.update({ "{}-{}".format(res.id[0],res.id[1]) : list(res['N'].coord)})
            except KeyError:
                resnum +=1
                continue
                
            try:
                O_coords.update({ "{}-{}".format(res.id[0],res.id[1]) : list(res['O'].coord)})
            except KeyError:
                resnum +=1
                continue
            resnum +=1
                
          
    return N_coords, O_coords



def get_vectors_ON(N_coords, O_coords):
    N_vector, O_vector = [[],[]]
    
    if N_coords != {}:
        N_vector, N_residues = protein_geometry.lstsq_determine_membrane_surface_plane(pdb_coords_per_gene_leaf=N_coords)
    
    
    if O_coords != {}:
        O_vector, O_residues = protein_geometry.lstsq_determine_membrane_surface_plane(pdb_coords_per_gene_leaf=O_coords)
    
    
    return N_vector, O_vector

def get_closest_leaf(N_dist, O_dist):
    if N_dist < O_dist:
        return 'N', 'O'
    else:
        return 'O', 'N'

def get_coords_of_membrane_crossing_residues(sfile, leaflet_residues):
    if type(sfile) == str:
        sfile = StructProp(ident = "temp",structure_path=sfile, file_type = op.basename(sfile).split('.')[-1])
        sfile = sfile.parse_structure()
    
    leaf_coords = {}
    for res in leaflet_residues:
        [chainId, resId] = res.split('_') 
        coord = sfile.first_model[chainId][int(resId)]['CA'].coord
        leaf_coords.update({res : list(coord)})
    return leaf_coords