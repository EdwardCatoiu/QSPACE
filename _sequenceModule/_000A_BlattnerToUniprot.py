# -*- coding: utf-8 -*-
from .utils import *
# import utils

log = logging.getLogger(__name__)



########################################################################################

def run_000A(geneListInput = '000A-gene_list_4349.txt',
        UniprotQueryInput = '000A-uniprot_ids_from_websearch.xlsx',
        manual_correction = False,
          trim = False,
       ):
    
    """
    ### Description:
    000A-SEQUENCE-blattner_to_uniprot_mapping

    The purpose of this notebook is to map gene IDs (blattner numbers) to uniprot IDs using the UniProt Search API
    (https://www.uniprot.org/id-mapping). Use "ensemble genomes" from UniprotKB

    ### Input(s)
    1. Gene list of Blattner IDs >> 000A-gene_list_4349.txt
    2. Search for the blattner IDs on UniProt and download the results as csv >> 000A-uniprot_ids_from_websearch.xlsx

    ### Output(s)
    1. {Blattner : UniProt ID} dictionary   >>  000A-uniprot_blattner_mapping.json 
    2. {Blattner : UniProt ID} dataframe    >>  000A-blattner_to_uniprot_mapping.csv

    ### Folder(s)
    -------------------------Inputs-------------------------
    qspace 
    └───inputs
    │   |   000A-gene_list_4349.txt 
    │   |   000A-uniprot_ids_from_websearch.xlsx
    |   |   ...
    -------------------------Outputs------------------------
    └───data
    │   |   000A-uniprot_blattner_mapping.json 
    │   |   000A-blattner_to_uniprot_mapping.csv
    |   |   ...


    """

    print ("Determining UniProt IDs from a list of Genes...\n-------------------------------------------------------------")
    infile = op.join(qspaceDirs["Input_dir"], geneListInput)
    

    f = open(infile)
    lines = f.readlines()
#     print lines[0:5]
    genelist = []
    for gene in lines:
        genelist +=[gene.split('\n')[0].split(' #')[0]]
    f.close()
#     print ("Input...")
    
    log.info('Gene list input file :\n\t{}'.format(infile))
 

    infile = op.join(qspaceDirs["Input_dir"], UniprotQueryInput)
    dfuniprot_websearch = pd.read_excel(infile, index_col = 0 )
    dfuniprot_websearch
    
    log.info('Uniprot query input file :\n\t{}'.format(infile))

#     print ("Calculating...")

    
    uniprot_mapping = dfuniprot_websearch['Entry'].to_dict()
 

    
    
    if manual_correction:
        uniprot_edit = {'b1370' : 'P76071','b2092' : 'P69831','b3643':'P0CG19','b2092':'P69831','b1370':'P76071',
                        'b3643':'P0CG19','b0218': 'P77354','b0229': 'Q47153', 'b0230': 'Q47154','b0235': 'P75675', 
                        'b0236': 'P28369', 'b0240': 'P24251', 'b0297': 'P36943', 'b0300': 'P77601', 'b0392': 'P75704',
                        'b0499': 'P77759', 'b0501': 'P33669', 'b0539': 'P75717', 'b0542': 'P75718', 'b0553': 'P21420',
                        'b0562': 'P77460', 'b0691': 'P37003', 'b0703': 'P77779', 'b0705': 'P75741', 'b1028': 'Q7DFV4',
                        'b1152': 'P75981', 'b1157': 'P33227','b1362': 'P77551', 'b1456': 'P24211',
                        'b1459': 'P76119', 'b1470': 'P76122', 'b1506': 'P76138', 'b1577': 'Q47138',
                        'b1578': 'P0CF60', 'b1579': 'P76168', 'b1936': 'P76323', 'b1999': 'P76359',
                        'b2139': 'P33369','b2190': 'P33924', 'b2238': 'P45505', 'b2355': 'P76508', 
                        'b2648': 'P76611', 'b2850': 'Q46786', 'b2854': 'Q46790', 'b2856': 'Q46791',
                        'b2858': 'Q46793', 'b2859': 'Q46795', 'b2862': 'Q46796', 'b2915': 'P64562', 
                        'b2969': 'Q46833',  'b2970': 'Q46834', 'b3046': 'P76655', 'b3134': 'P42905', 
                        'b3135': 'P42906', 'b3268': 'P45766','b3423': 'P0ACL0', 'b3443': 'P46856',
                        'b3504': 'P37635',  'b3534': 'P37655', 'b3595': 'P32109',  'b4038': 'P32690',
                        'b4205': 'P39309', 'b4271': 'P39347', 'b4281': 'P39354', 'b4282': 'P39355',
                        'b4286': 'Q47719', 'b4308': 'P39369', 'b4462': 'P76616', 'b4486': 'P39393',
                        'b4488': 'P0DP89',  'b4495': 'P76322', 'b4497': 'P76349','b4503': 'Q2EEP9',
                        'b4521': 'P76000', 'b4543': 'V9HVX0', 'b4569': 'P45421', 'b4570': 'P77184',
                       'b4573': 'P30192',  'b4580': 'P77199', 'b4581': 'P0DP63',  'b4583': 'P0DP69',
                        'b4584': 'P0DP21','b4587': 'P75679', 'b4623': 'P39212',  'b4658': 'Q46850',
                        'b4660': 'P37629',  'b4696': 'P32051','b0115' : 'P06959','b1474':'P24183',
                        'b4079':'P07658','b4403':'P37005','b3894':'P32176', 
                        
                       }

        uniprot_mapping.update(uniprot_edit)   #EDIT,
        log.info('Using manual UniprotID mapping'.format(len(uniprot_edit)))

        
    if trim:
        dellist = []
        for k in uniprot_mapping:
            if k not in genelist:
                dellist +=[k]
        for k in dellist:
            del uniprot_mapping[k]

        
    log.info('Found UniProt IDs for {} genes'.format(len(uniprot_mapping)))
#     print ("Output...")

    outfile = op.join(qspaceDirs["DataOutput_dir"],'000A-uniprot_blattner_mapping.json')
#     print ("\t> {}".format(outfile))
    log.info('Saving Blattner-Uniprot mapping :\n\t{}'.format(outfile))

    with open(outfile, 'w') as f:
        json.dump(uniprot_mapping, f)

    df = pd.DataFrame.from_dict(uniprot_mapping, orient = 'index')
    outfile = op.join(qspaceDirs["DataOutput_dir"],'000A-blattner_to_uniprot_mapping.csv')
    log.info('Saving Blattner-Uniprot mapping :\n\t{}'.format(outfile))

    df.to_csv(outfile)
    return uniprot_mapping, genelist