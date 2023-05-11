# -*- coding: utf-8 -*-

from .utils import *
log = logging.getLogger(__name__)


def run_001D(genelist, blattnerUniprotMapping ,alphafold = True, swiss = True, itasser = True,itasser_metadata=False,force_rerun = False):
    """
    Sequence Alignment of AAseqs to Homology models 
    
    """
    
    def getDfSeq(genelist=genelist, blattner_uniprot_mapping=blattnerUniprotMapping):
    
    
        dfseq= pd.DataFrame()
        print ("Getting gene sequences from Uniprot and Alleleome...\n---------------------------------------------")
        for blattner in tqdm(set(genelist).union(set(blattner_uniprot_mapping.keys()))):
            if blattner in blattner_uniprot_mapping:
                uniprotId = blattner_uniprot_mapping[blattner]
                uniprotFasta = op.join(qspaceDirs['UniprotSeqsDir'], '{}.fasta'.format(uniprotId))
                dfseq.loc[blattner, 'UniProtId'] = uniprotId
                if op.exists(uniprotFasta):
                    record_uni = SeqIO.read(uniprotFasta, "fasta")
                    dfseq.loc[blattner, 'UniProtSeq'] = str(record_uni.seq)
                    dfseq.loc[blattner, 'UniProtSeqLen'] = len(str(record_uni.seq))



            alleleomeFasta =op.join(qspaceDirs['AlleleomeSeqsDir'],'{}_WT_consensus.fasta'.format(blattner))       
            if op.exists(alleleomeFasta):
                record_wt = SeqIO.read(uniprotFasta, "fasta")
                dfseq.loc[blattner, 'AlleleomeId'] = '{}_WT_consensus'.format(blattner)
                dfseq.loc[blattner, 'AlleleomeSeq'] = str(record_wt.seq)
                dfseq.loc[blattner, 'AlleleomeSeqLen'] = len(str(record_wt.seq))
        return dfseq
    
    dfseq = getDfSeq(genelist=genelist, blattner_uniprot_mapping=blattnerUniprotMapping)
    outfile = op.join(qspaceDirs['DataOutput_dir'], '001D-dfrepseq.csv')
    dfseq.to_csv(outfile)
    log.info('Saving DataFrame of gene sequences ... \n\t{}'.format(outfile))

#     print "Getting representative sequences..."
#     print "\t> {}".format(outfile)
    
    
    if swiss:
        infile  = op.join(qspaceDirs['DataOutput_dir'], '001B-expected_alignment_files_SWISS.json')
        with open(infile, 'r') as f:
            expected_files = json.load(f)
        
        #total Chains
        infile  = op.join(qspaceDirs['DataOutput_dir'], '001B-SWISS-model_chains_to_genes.json')
        with open(infile, 'r') as f:
            total_chains = json.load(f)
    
   
        alignment_outfolder = qspaceDirs['SequenceAlignmentDir']
        checked = []
        swiss_structures_folder =  qspaceDirs['swissStructuresDir']
        cif_folder = qspaceDirs['cifStructuresDir']
        error_needle = []
        structure_error = []

        print ('Aligning SWISS-PROT models to Uniprot Sequences...\n------------------------------------------------------------')
        for pdb_id , mfiles in tqdm(expected_files.items()):
        #     print pdb_id

            #remove consensus sequences for now
            missing_files = []
            for m in mfiles:
                if 'consensus' in m:
                    continue
                missing_files += [m]

        #     print uniprot_missing_files
            checked ,structure_error,error_needle = alignment.align_missing_files(pdb_id,                                                   
                                                                                  missing_needle_files=missing_files,
                                                                                  dfrepseq=dfseq,
                                                                                  structure_error =structure_error ,
                                                                                  error_needle= error_needle,
                                                                                  force_rerun= force_rerun,
                                                                                  struch_seq_alignment_folder=alignment_outfolder,
                                                                                  checked_pdbs= checked, 
                                                                                  pdb_folder=swiss_structures_folder,
                                                                                  cif_folder=cif_folder,
                                                                                  using_swiss = True
                                                                                 )
            
        input_quality = {}
        c = 0 
        m = 0 
#         print 'Getting SWISS models sequence identity...'

        for pdb_id, gene_maps in tqdm(total_chains.items()):

            temp = {}
            for gene, chain_list in gene_maps.items():
                info = alignment.get_input_quality(pdb_id = pdb_id,
                                                   gene= gene,
                                                   chain_list =chain_list, 
                                                   dfrepseq= dfseq,
                                                   use_rep = True,
                                                   alignment_folder= qspaceDirs['SequenceAlignmentDir'])
                if info == {gene : {}}:
                    c +=1
#                     print pdb_id, gene,'\t', c,'\t', m
                    continue

                temp.update(info)
                m +=1
            input_quality.update({pdb_id : temp})
        outfile = op.join(qspaceDirs['DataOutput_dir'], '001D-quality_of_swiss_models.json')
        with open(outfile, 'w') as f:
            json.dump(input_quality,f)
        log.info("Saving SWISS-PROT model quality ... \n\t{}".format(outfile))


            
    
    if itasser:

        #missing files
        infile  = op.join(qspaceDirs['DataOutput_dir'], '001C-expected_alignment_files_ITASSER.json')
        with open(infile, 'r') as f:
            expected_files = json.load(f)
    
        #total Chains
        infile  = op.join(qspaceDirs['DataOutput_dir'], '001C-ITASSER-model_chains_to_genes.json')
        with open(infile, 'r') as f:
            total_chains = json.load(f)
            
   
        alignment_outfolder = qspaceDirs['SequenceAlignmentDir']
        checked = []
        structures_folder =  qspaceDirs['itasserCleanStringRemovedStructuresDir']
        cif_folder = qspaceDirs['cifStructuresDir']
        error_needle = []
        structure_error = []

        print ('Aligning I-TASSER models to Uniprot Sequences...\n------------------------------------------------------------')
        for pdb_id , mfiles in tqdm(expected_files.items()):
            pdb_id = op.basename(pdb_id)
        #     print pdb_id

            #remove consensus sequences for now
            missing_files = []
            for m in mfiles:
                if 'consensus' in m:
                    continue
                missing_files += [op.basename(m)]

        #     print missing_files
            checked ,structure_error,error_needle = alignment.align_missing_files(pdb_id, 
                                                                                  missing_needle_files=missing_files,
                                                                                  dfrepseq=dfseq,
                                                                                  structure_error =structure_error ,
                                                                                  error_needle= error_needle,
                                                                                  force_rerun= force_rerun,
                                                                                  struch_seq_alignment_folder=alignment_outfolder,
                                                                                  checked_pdbs= checked, 
                                                                                  pdb_folder=structures_folder,
                                                                                  cif_folder=cif_folder,
                                                                                  using_swiss = False,
                                                                                  using_itasser = True
                                                                                 )
            
            
            
#         print 'Getting ITASSER models sequence identity...'
        input_quality = {}
        c = 0 
        m = 0 
        for pdb_id, gene_maps in tqdm(total_chains.items()):

            temp = {}
            for gene, chain_list in gene_maps.items():
                info = alignment.get_input_quality(pdb_id = pdb_id,
                                                   gene= gene,
                                                   chain_list =chain_list, 
                                                   dfrepseq= dfseq,
                                                   use_rep = True,
                                                   alignment_folder= qspaceDirs['SequenceAlignmentDir'])
                if info == {gene : {}}:
                    c +=1
#                     print pdb_id, gene,'\t', c,'\t', m
                    continue

                temp.update(info)
                m +=1
            input_quality.update({pdb_id : temp})
            
            
        input_quality_tm_score = {}
#         print 'Getting ITASSER models quality (sequence identity x tm-score)...'

        for itasserId, input_quality in tqdm(input_quality.items()):
            itasserEntry = itasserId.split('_model')[0]
            tm_Score = itasser_metadata[itasser_metadata.get('Entry Name') == itasserEntry].get('TM-score').values[0]


            for gene, chainInfo in input_quality.items():
                for chain, chainQualInfo in chainInfo.items():
                    [bool_metric, qual_metric] = chainQualInfo
                    new_qual_metric = qual_metric * tm_Score

                    if itasserId not in input_quality_tm_score:
                        input_quality_tm_score.update({itasserId : {gene : {chain : [bool(new_qual_metric > 50), new_qual_metric]}}})

                    elif gene not in input_quality_tm_score[itasserId]:
                        input_quality_tm_score[itasserId].update({gene : {chain : [bool(new_qual_metric > 50), new_qual_metric]}})
                    elif chain not in input_quality_tm_score[itasserId][gene]:
                        input_quality_tm_score[itasserId][gene].update({chain : [bool(new_qual_metric > 50), new_qual_metric]})
        outfile = op.join(qspaceDirs['DataOutput_dir'], '001D-quality_of_itasser_models.json')
        with open(outfile, 'w') as f:
            json.dump(input_quality_tm_score,f)
        log.info("Saving ITASSER model quality ... \n\t{}".format(outfile))
    
    
    if alphafold:

        infile  = op.join(qspaceDirs['DataOutput_dir'], '001A-expected_alignment_files_ALPHAFOLD.json')
        with open(infile, 'r') as f:
            expected_files = json.load(f)
            
        
    
    
   
        alignment_outfolder = qspaceDirs['SequenceAlignmentDir']
        checked_itasser = []
        structures_folder =  qspaceDirs['alphafoldStructuresDir']
        cif_folder = qspaceDirs['cifStructuresDir']
        error_needle = []
        structure_error= []

        print ('Aligning AlphaFold models to Uniprot Sequences...\n------------------------------------------------------------')
        for pdb_id , mfiles in tqdm(expected_files.items()):
            pdb_id = op.basename(pdb_id)
        #     print pdb_id

            #remove consensus sequences for now
            missing_files = []
            for m in mfiles:
                if 'consensus' in m:
                    continue
                missing_files += [op.basename(m)]

        #     print missing_files
            checked ,structure_error,error_needle = alignment.align_missing_files(pdb_id,
                                                                                  missing_needle_files=missing_files,
                                                                                  dfrepseq=dfseq,
                                                                                  structure_error =structure_error ,
                                                                                  error_needle= error_needle,
                                                                                  force_rerun= force_rerun,
                                                                                  struch_seq_alignment_folder=alignment_outfolder,
                                                                                  checked_pdbs= checked, 
                                                                                  pdb_folder=structures_folder,
                                                                                  cif_folder=cif_folder,
                                                                                  using_swiss = False,
                                                                                  using_itasser = False,
                                                                                  using_alphafold = True
                                                                                 )

    return dfseq
