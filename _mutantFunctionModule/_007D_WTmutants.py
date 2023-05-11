# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)


def run_007D_mapWTmutants(dfalldata,
                          alleleomeInputFolder='/home/ecatoiu/Projects/Nature_Alleleome/Alleleome_Data/dfz',
                          force_realign = True,
                          
                         ):

    global_index_dfz_index = {}
    for cplx in tqdm(dfalldata.cplx.unique()):
        dfcplx = dfalldata[dfalldata.cplx == cplx]
        for structureId in (dfcplx.structureId.unique()):
            dfs = dfcplx[dfcplx.structureId == structureId]
            for chain in dfs.sfileChainId.unique():
                dfchain = dfs[dfs.sfileChainId ==chain]
                for gene in dfchain.gene.unique():
                    dfg = dfchain[dfchain.gene == gene]
                    dfg = dfg[dfg.seqAA1.isna() == False]
                    glob_index_list = dfg.index.tolist()

                    infile = op.join(alleleomeInputFolder, '{}_dfz.csv'.format(gene))
                    try:
                        dfz = pd.read_csv(infile, index_col=0)
                        
                    except IOError:
                        log.info("No WT variation file. Affects : {} {} {}".format(cplx, structureId , infile))
                        continue
                    dfz=dfz[dfz.consensus_AAseq.isna() == False]
                    dfz=dfz[dfz.consensus_AAseq.isin(['-','*']) == False]
                    dfz_index_list  = dfz.index.tolist()

                    seq_structproteome = "".join(dfg.seqAA1.tolist())
                    seq_alleleome = "".join(dfz.consensus_AAseq.tolist())

                    needle_file = ssbioaln.run_needle_alignment(seq_a = seq_structproteome,seq_b = seq_alleleome,force_rerun=force_realign)
                    dfaln = ssbioaln.get_alignment_df_from_file(needle_file)
                    temp = dfaln[dfaln.id_a_aa.isna() == False]
                    temp['global_index_a'] = pd.Series(glob_index_list)
                    dfaln['global_index_a'] = temp['global_index_a']

                    temp = dfaln[dfaln.id_b_aa.isna() == False]
                    temp['dfz_index_b'] = pd.Series(dfz_index_list)
                    dfaln['dfz_index_b'] = temp['dfz_index_b']

                    dfaln_glob_index = dfaln.global_index_a.to_dict()
                    dfaln_dfz_index = dfaln.dfz_index_b.to_dict()

                    indexmap = {dfaln_glob_index[i] : dfaln_dfz_index[i] for i in dfaln_glob_index.keys() }
                    global_index_dfz_index.update(indexmap)
    outfile = op.join(qspaceDirs['DataOutput_dir'], '007D-StructuralProteome_global_dfz_index_WT_variation.json') 
    with open(outfile,'w') as f:   
        json.dump(global_index_dfz_index, f)
    print ("Saving....\n\t> {}".format(outfile))
    
    
    global_gene_index = dfalldata.gene.to_dict()
    glob_WT_AAseq_details = {}
    glob_WT_AAseq_dominant_freq = {}
    glob_WT_CODONseq_details = {}
    glob_WT_CODONseq_dominant_freq = {}

    for g in tqdm(dfalldata.gene.unique()):
        dfg = dfalldata[dfalldata.gene == g]
#         print g
        dfz = op.join(alleleomeInputFolder, '{}_dfz.csv'.format(g))
        try:
#             print g, dfz
            dfz = pd.read_csv(dfz, index_col=0)
            
        except IOError:
            print (cplx, structureId, chain, gene , '<<<DNE')
            continue

        WT_AAseq_details = dfz.consensus_AAseq_details.to_dict()
        WT_AAseq_dominant_freq = dfz.consensus_strains_normalized.to_dict()
        WT_CODONseq_details = dfz.consensus_codon_seq_details.to_dict()
        WT_CODONseq_dominant_freq = dfz.consensus_strains_normalized_codon.to_dict()

        for index in dfg.index.tolist():
    #         print index
            if index in global_index_dfz_index:
                dfz_index = global_index_dfz_index[index]
                if str(dfz_index) == 'nan':
                    continue
                glob_WT_AAseq_details.update({index : WT_AAseq_details[dfz_index]})
                glob_WT_AAseq_dominant_freq.update({index : WT_AAseq_dominant_freq[dfz_index]})
                glob_WT_CODONseq_details.update({index : WT_CODONseq_details[dfz_index]})
                glob_WT_CODONseq_dominant_freq.update({index : WT_CODONseq_dominant_freq[dfz_index]})

                
    outfile = op.join(qspaceDirs['DataOutput_dir'], '007D-StructuralProteome_WT_AA_seq_details.json') 
    with open(outfile,'w') as f:   
        json.dump(glob_WT_AAseq_details, f)
    print ("Saving....\n\t> {}".format(outfile))
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '007D-StructuralProteome_WT_AA_dom_frequency.json') 
    with open(outfile,'w') as f:   
        json.dump(glob_WT_AAseq_dominant_freq, f)
    print ("Saving....\n\t> {}".format(outfile))
    
    outfile = op.join(qspaceDirs['DataOutput_dir'], '007D-StructuralProteome_WT_CODON_seq_details.json') 
    with open(outfile,'w') as f:   
        json.dump(glob_WT_CODONseq_details, f)
    print ("Saving....\n\t> {}".format(outfile))
                
    outfile = op.join(qspaceDirs['DataOutput_dir'], '007D-StructuralProteome_WT_CODON_dom_frequency.json') 
    with open(outfile,'w') as f:   
        json.dump(glob_WT_CODONseq_dominant_freq, f)
    print ("Saving....\n\t> {}".format(outfile))
                
    return global_index_dfz_index, glob_WT_CODONseq_details, glob_WT_CODONseq_dominant_freq, glob_WT_AAseq_details, glob_WT_AAseq_dominant_freq
