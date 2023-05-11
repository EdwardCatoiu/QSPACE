import os.path as op
import copy
import os
import pandas as pd
import sys
import check_structure

from ssbio.protein.structure.structprop import StructProp
import ssbio.protein.sequence.utils.alignment as ssbioaln




def find_expected_needle_files(pdb_id, 
                               total_chains = {},
                               blattner_uniprot_mapping = {},
                               structure_error=[],
                               error_needle = [],
                               force_rerun = False,
                               seq_alignment_folder = '../GEMPRO/needle_alignments/'
                              ):
#long cut (DOWNLOAD ALL THE ALLIGNMENTS)
# pdbs_to_check = list(set(total_chains.keys()) - set(error_list) - set(checked_pdbs))

    if not seq_alignment_folder:
        seq_alignment_folder = '../GEMPRO/needle_alignments/'

    
    genelist = total_chains[pdb_id].keys()
    expected_needle_files = []
    for gene in genelist:
        repseq_id = blattner_uniprot_mapping[gene]
        consensus_id = 'consensus_allele'       

        mapped_chains = total_chains[pdb_id][gene]
        for chain_id in mapped_chains:
            aln_id1 = '{}_{}_{}-{}'.format(pdb_id, gene, repseq_id, chain_id)
            outfile1 = "{}.needle".format(aln_id1)
            aln_id2 = '{}_{}_{}-{}'.format(pdb_id, gene, consensus_id, chain_id)
            outfile2 = "{}.needle".format(aln_id2)
            
            expected_needle_files +=[outfile1]
            expected_needle_files +=[outfile2]
        
    missing_needle_files = []
    for needle_file in expected_needle_files:
        if not op.exists(op.join(seq_alignment_folder, needle_file)):
            missing_needle_files +=[needle_file]
    missing_needle_files = list(set(missing_needle_files) - set(error_needle))
    
    if len(missing_needle_files) == 0 and not force_rerun:
#         print 'FINE : {}'.format(pdb_id)
        return missing_needle_files
    
    elif force_rerun:
        missing_needle_files = copy.deepcopy(expected_needle_files)
        return missing_needle_files
   
    return missing_needle_files


def align_missing_files(pdb_id, 
                        missing_needle_files,
                        pdb_folder,
                        cif_folder,
                        struch_seq_alignment_folder,
                        dfrepseq,
                       structure_error = [],
                        error_needle= [],
                        force_rerun= False,
                        checked_pdbs  = [],
                        using_swiss = False,
                        using_itasser = False,
                        using_alphafold = False
                       ):
    
    if not force_rerun:
        missing_needle_files = list(set(missing_needle_files) - set(error_needle))
    
    if len(missing_needle_files) ==0:
        checked_pdbs +=[pdb_id]
        return checked_pdbs, structure_error,error_needle
    
    if pdb_id in checked_pdbs + structure_error:
        return checked_pdbs,structure_error,error_needle
    
    if not using_swiss and not using_itasser and not using_alphafold:
        pdb_id = pdb_id.lower()
    
    #load the structure file
    cif_file = op.join(cif_folder, "{}.cif".format(pdb_id))
    pdb_file = op.join(pdb_folder, "{}.pdb".format(pdb_id))
#     print op.exists(pdb_file), pdb_file
    if op.exists(pdb_file):
        #print 'pdb file exists'
        a_structure = StructProp(ident = pdb_id,structure_path=pdb_file,file_type = 'pdb')
#         print pdb_file
    elif op.exists(cif_file):
        a_structure = StructProp(ident = pdb_id,structure_path=cif_file,file_type = 'cif')
#         print cif_file
    try:
        parsed = a_structure.parse_structure()
    except (AttributeError, IndexError, PDBConstructionException):
#         log.info( "{} Structure Error".format(pdb_id))
        structure_error +=[pdb_id]
        return    checked_pdbs, structure_error,error_needle
    
    #print pdb_id,
    for needle_file in missing_needle_files:
        if needle_file in error_needle:
            continue
#         aln_id = '{}_{}_{}-{}'.format(pdb_id, gene, repseq_id, chain_id)
        
        if not using_swiss and not using_itasser and not using_alphafold:
            pdb_id = needle_file.split('_')[0]
            gene = needle_file.split('_')[1]
            chain_id = needle_file.split('.needle')[0].split('-')[1]
#             repseq_id = needle_file.split('{}_'.format(gene))[1].split('-')[0]
        elif using_swiss:
            chain_id = needle_file.split('.needle')[0].split('-')[1]
            gene = needle_file.split('_')[4]
#             repseq_id = needle_file.split('{}_'.format(gene))[1].split('-')[0]
            pdb_id = needle_file.split('_{}'.format(gene))[0]
        elif using_itasser:
            if 'ECOLI' not in needle_file:
                chain_id = needle_file.split('.needle')[0].split('-')[1]
                gene = needle_file.split('_')[4]
#                 repseq_id = needle_file.split('{}_'.format(gene))[1].split('-')[0]
                pdb_id = needle_file.split('_{}'.format(gene))[0]

            else:
                chain_id = needle_file.split('.needle')[0].split('-')[1]
                gene = needle_file.split('_')[6]
#                 repseq_id = needle_file.split('{}_'.format(gene))[1].split('-')[0]
                pdb_id = needle_file.split('_{}'.format(gene))[0]
        elif using_alphafold:
            chain_id = needle_file.split('.needle')[0].split('-')[-1]   
            gene = needle_file.split('{}_'.format(pdb_id))[-1].split('_')[0]


            
        
        
        if 'consensus' in needle_file:
            g_seq = dfrepseq.loc[gene,'AlleleomeSeq']
        else:
            g_seq = dfrepseq.loc[gene,'UniProtSeq']
        
        if g_seq is None or str(g_seq) == 'nan' :
#             print ( '\t Error : No Sequence >>>> {}' .format(needle_file))
            error_needle +=[needle_file]
            continue
               
        chain_seq_record = None
#         print 'GOT HERE'
#         print chain_id
        try:
             chain_prop = a_structure.chains.get_by_id(chain_id)
             chain_seq_record = chain_prop.seq_record
        except (KeyError,AttributeError):
#             print ( '{}-{} Chain Does Not Exist'.format(pdb_id, chain_id))
            error_needle +=[needle_file]
            continue
        if not chain_seq_record:
#             print  ( "Chain Sequence Not Parsed : {}-{}".format(pdb_id, chain_id))
            error_needle +=[needle_file]
            continue 
        if len(chain_seq_record) == 0:
#             print ( "Not a real chain :  {}-{}".format(pdb_id, chain_id))
            error_needle +=[needle_file]
            continue       
        
       
        #SEQUENCES 
        struct_seq = str(chain_prop.seq_record.seq)
        if struct_seq.count('X') > 0:
            #eliminate the HOH X in the chain seq
            struct_seq = "".join([''.join(aa) for aa in str(chain_prop.seq_record.seq) if aa != 'X'])
            if len(struct_seq) == 0:
                error_needle +=[needle_file]
                continue
#         print '\t Test >>>> {}'.format( needle_file)

        outfile_path = op.join(struch_seq_alignment_folder, needle_file)
        if not op.exists(outfile_path) or force_rerun:
            if type(g_seq) == str:
                aln_file = ssbioaln.run_needle_alignment(seq_a = struct_seq, seq_b = g_seq,
                                                         outfile = needle_file, 
                                                         outdir = struch_seq_alignment_folder,
                                                         force_rerun=force_rerun)
#                 log.info('>>>> {}'.format( needle_file))
    checked_pdbs +=[pdb_id]
    return    checked_pdbs, structure_error,error_needle

def get_input_quality(pdb_id, gene, chain_list, dfrepseq, alignment_folder ,use_rep = True):
    
    if use_rep:
        repseq_id = dfrepseq.loc[gene,'UniProtId']
    else:
        repseq_id = 'consensus_allele'
    temp = {}
    for chain_id in chain_list:
        needle_file = '{}_{}_{}-{}.needle'.format(pdb_id, gene, repseq_id, chain_id)    
        needle_path = op.join(alignment_folder, needle_file)
        if not op.exists(needle_path):
            continue
        needle_stats = ssbioaln.needle_statistics(needle_path)['asis_asis']
        if needle_stats =={}:
#             print needle_file, 'no stats'
            continue
        percent_identity = needle_stats['percent_identity']
        info = [percent_identity >= 50, percent_identity]
        temp.update({chain_id : info})
    
    return {gene : temp}

