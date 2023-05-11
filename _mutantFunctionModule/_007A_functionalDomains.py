# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)


def run_007A_functionalFeatures(dfalldata,
                                force_rerun = True,
                                uniprotSeqDir = qspaceDirs['UniprotSeqsDir'],
                                featuresUniprotDir = qspaceDirs['featuresUniprotDir'],
                               ):
    print ("Parsing Functional Annotations from UniProt ...\n------------------------------------")

    count = 0 
    uniprotIds = dfalldata[dfalldata.geneSeqId !='WT_consensus'].geneSeqId.unique()
    
    for uniprotId in tqdm(uniprotIds):
        outFileCsv = op.join(featuresUniprotDir,'{}_features.csv'.format(uniprotId))
        if op.exists(outFileCsv) and not force_rerun:
            count +=1
            continue 

        textpath = op.join(uniprotSeqDir, '{}.txt'.format(uniprotId))
        fastapath = op.join(uniprotSeqDir, '{}.fasta'.format(uniprotId))

        feature_df =  uniprot_features.make_nextGenome_feature_dataframe(uniprotId,
                                                                         uniprot_text_path=textpath, 
                                                                         uniprot_fasta_path=fastapath
                                                                        )
        if len(feature_df) != 0:
            feature_df.to_csv(outFileCsv)
            count += 1
            
            
    print ("Mapping Functional Annotations to the structural proteome ...\n------------------------------------")
      
    global_feature_dict = {}

    for uniprot_df_file in tqdm(os.listdir(featuresUniprotDir)):
    #     print uniprot_df_file
        f = pd.read_csv(op.join(featuresUniprotDir, uniprot_df_file),index_col = 0 )

        uniprot_id = uniprot_df_file.split('_')[0]
        if uniprot_id not in dfalldata.geneSeqId:
            continue
        
        t = dfalldata[dfalldata.geneSeqId  == uniprot_id]
        
        for feat, row in f.iterrows():
            if feat in ['AA_Seq', 'AA_SeqNum','Gene_ID','UniProt']:
                continue
    #         if feat in ['MUTAGEN' , 'CONFLICT','VARIANT']:
    #             continue
    #         if feat in ['HELIX', 'TURN','STRAND','TRANSMEM', 'TOPO_DOM']:
    #             continue

            if feat not in global_feature_dict:
                global_feature_dict.update({feat : {}})

            region = f[f.index == feat].dropna(axis = 1)
            for c in region.columns:
                feature_info = region.loc[feat, c]

                seq_resnum = int(c.split('_')[1])
    #             if feature_info =="['Unknown']":
    #                 print uniprot_id, seq_resnum, feat
                for index_all in t[t.seqNum == seq_resnum].index.tolist():


                    global_feature_dict[feat].update({index_all : feature_info})
                
    outfile = op.join(qspaceDirs['DataOutput_dir'],'007A-StructuralProteome_UniprotFeatures_dict.json')    
    with open(outfile, 'w') as f:
        json.dump(global_feature_dict,f )
    log.info("Saving QSPACE Proteome - UniProt Functional Annotations...\n\t{}".format(outfile))
    return global_feature_dict            
                        