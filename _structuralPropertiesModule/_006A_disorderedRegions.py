# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)


def disordered_disembl_regions(SeqRecord,disembl_cmd):
    
    def zeroorone(x):
        import string
        if x in string.ascii_uppercase:
            return 1
        else:
            return 0

    with tempfile.NamedTemporaryFile(delete=True) as f:
        SeqIO.write(SeqRecord, f.name, "fasta")
        out = subprocess.check_output(
            ['python2', disembl_cmd, '8', '8', '4', '1.2', '1.4', '1.2', f.name]
        )
        out2 = out.decode().split('\n')
        coils = list(map(lambda x: zeroorone(x), out2[9]))
        rem465 = list(map(lambda x: zeroorone(x), out2[11]))
        hotloops = list(map(lambda x: zeroorone(x), out2[13]))

        
        return coils,rem465,hotloops
    
    
def run_006A_Disembl(dfalldata , disembl_cmd,dfFastaLoc = pd.DataFrame()):
    
    def parse_DISEMBL(gene_id , df_disembl,df):
    
        df_disemble_gene = df_disembl[df_disembl.index == gene_id]
        disemble_dict = {'disembl_coils':{},'disembl_rem465':{},'disembl_hotloops' :{}}


        for feat in disemble_dict:
            feature_info_list = ast.literal_eval(df_disemble_gene[feat].values[0])
            info = map(lambda index, feat_info : {int(index) : feat_info}, df.index.values, feature_info_list)
            for i in info:
                disemble_dict[feat].update(i)
        return disemble_dict

    print ("Running DisEMBL...")
    if len(dfFastaLoc) ==0:
           dfFastaLoc = get_fastaFileLocations(dfalldata,
                                              alleleomeDir = qspaceDirs['AlleleomeSeqsDir'],
                                              uniprotDir = qspaceDirs['UniprotSeqsDir'],
                                             )
    
    df_disemble = copy.deepcopy(dfFastaLoc)
    index_list = df_disemble.index.tolist()

    for index in tqdm(index_list):
        row = df_disemble.loc[index]
        seqfile = df_disemble.loc[index]['fasta_location']
        fasta_id = df_disemble.loc[index]["fasta_id"]
#         print (seqfile)
        record = SeqIO.read(seqfile, "fasta")

#         seqprop = SeqProp(id =fasta_id,seq = str(record.seq), sequence_path =seqfile)

#         seqprop = SeqProp(ident =fasta_id, sequence_file=seqfile)
        coils,rem465,hotloops = disordered_disembl_regions(SeqRecord=record,disembl_cmd=disembl_cmd)
        df_disemble.loc[index, 'disembl_coils'] = str(coils)
        df_disemble.loc[index, 'disembl_rem465'] = str(rem465)
        df_disemble.loc[index, 'disembl_hotloops'] = str(hotloops)
    outfile = op.join(qspaceDirs['DataOutput_dir'],'006A_DISEMBL_dataframe.csv')
    df_disemble.to_csv(outfile)
    log.info("Saving DisEMBL Proteins dataframe...\n\t{}".format(outfile))
    
    print ('AA-level disorder annotation...')
    cplx_list = dfalldata.cplx.unique()
    disemble_global_info = {'disembl_coils':{},'disembl_rem465':{},'disembl_hotloops' :{}}
    for cplx in tqdm(cplx_list):
        dfc = dfalldata[dfalldata.cplx ==cplx]
        for structure_id in dfc.structureId.unique():
            dfs = dfc[dfc.structureId ==structure_id]
            for gene_id in dfs.gene.unique():
                dfg = dfs[dfs.gene ==gene_id]
                for chain_id in dfg.get('structChainId_mod').unique():
                    dfgc = dfg[dfg.get('structChainId_mod') ==chain_id]
                    disemble_info = parse_DISEMBL(gene_id, df_disemble, df = dfgc)
                    for k , v in disemble_info.items():
                        disemble_global_info[k].update(v)
   
    outfile = op.join(qspaceDirs['DataOutput_dir'],'006A-StructuralProteome_DISEMBL_disorder.json')    
    with open(outfile, 'w') as f:
        json.dump(disemble_global_info,f )
    log.info("Saving QSPACE Proteome - DisEMBL AA annotations...\n\t{}".format(outfile))

    return df_disemble,  disemble_global_info