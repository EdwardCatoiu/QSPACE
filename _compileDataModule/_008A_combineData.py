# -*- coding: utf-8 -*-
from .utils import *
log = logging.getLogger(__name__)

def convert_index_to_int(data):
    data_int = {}
    for data_key,  dataInfo in tqdm(data.items()):
        for index, value in dataInfo.items():
            try:
                int(index)
            except ValueError:
                continue
            
            if data_key not in data_int:
                data_int.update({data_key : {int(index) : value}})
            elif int(index) not in data_int[data_key]:
                data_int[data_key].update({int(index) : value})
    return data_int

def run_008A_finalizeDataframe(dfalldata):
    
    #Protein Compartments
    filename = '005F-StructuralProteome_Protein_compartment.json'
    filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
    with open(filepath, 'r') as f:
        data = json.load(f)
    data_int = {}
    for k ,v in data.items():
        data_int.update({int(k) : v})
    dfalldata['005F_Protein_Compartment'] = pd.Series(data_int)
    log.info("Added .. {}".format(filename))
    
    #AA Compartments
    filename = '005F-StructuralProteome_AA_compartment.json'
    filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
    with open(filepath, 'r') as f:
        data = json.load(f)
    data_int = {}
    for k ,v in data.items():
        data_int.update({int(k) : v})
    dfalldata['005F_AA_Compartment'] = pd.Series(data_int) 
    log.info("Added .. {}".format(filename))

    #Protein Disordered Regions
    filename = '006A-StructuralProteome_DISEMBL_disorder.json'
    filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
    with open(filepath, 'r') as f:
        data = json.load(f)
    data_int = convert_index_to_int(data)
    for key, series in data_int.items():
        dfalldata["006A_{}".format(key)] = pd.Series(series)
    log.info("Added .. {}".format(filename))

        
    #SCRATCH
    filename =   '006B-StructuralProteome_SCRATCH.json'
    filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
    with open(filepath, 'r') as f:
        data = json.load(f)
    data_int = convert_index_to_int(data)
    for key, series in data_int.items():
        dfalldata["006B_SCRATCH_{}".format(key)] = pd.Series(series)
    log.info("Added .. {}".format(filename))
        
    #MSMS
#     filename =    '006C-StructuralProteome_MSMS_res_depth.json'
#     with open('../data/raw/{}'.format(filename), 'r') as f:
#         data = json.load(f)
#     data_int = {}
#     for k ,v in data.items():
#         if str(k) in ['nan','NaN','null']:
#             continue
#         data_int.update({int(float(k)) : v})
#     dfalldata["006C_MSMS_residue_depth"] = pd.Series(data_int)
    
#     filename =    '006C-StructuralProteome_MSMS_ca_depth.json'
#     with open('../data/raw/{}'.format(filename), 'r') as f:
#         data = json.load(f)
#     data_int = {}
#     for k ,v in data.items():
#         if str(k) in ['nan','NaN','null']:
#             continue
#         data_int.update({int(float(k)) : v})
#     dfalldata["006C_MSMS_alphacarbon_depth"] = pd.Series(data_int)

    #Disulfide Bridges
    filename =   '006D-StructuralProteome-disulfide_bridges.json'
    filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
    with open(filepath, 'r') as f:
        data = json.load(f)
    data_int = {}
    for k ,v in data.items():
        if str(k) in ['nan','NaN']:
            continue
        data_int.update({int(float(k)) : v})
    dfalldata["006D_Disfulfide_Bridges_3A"] = pd.Series(data_int)
    log.info("Added .. {}".format(filename))

    #protein protein interfaces (ScanNet)
    filename =  '006E-StructuralProteome-interfaces_ScanNet.json'
    filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
    with open(filepath, 'r') as f:
        data = json.load(f)
    data_int = {}
    for k ,v in data.items():
        if str(k) in ['nan','NaN']:
            continue
        data_int.update({int(float(k)) : v})
    dfalldata["006E_ppInterface_ScanNet"] = pd.Series(data_int)
    log.info("Added .. {}".format(filename))

        
    #protein protein interfaces (geometric)
    filename =  '006E-StructuralProteome-interfaces_geometry.json'
    filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
    with open(filepath, 'r') as f:  
        data = json.load(f)
    data_int = convert_index_to_int(data)
    for key, series in data_int.items():
        dfalldata["006E_ppInterface_Geometry_{}Angs".format(key)] = pd.Series(series)
    log.info("Added .. {}".format(filename))

        
    #DSSP 
    filename =  '006F-StructuralProteome_DSSP.json'
    filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
    with open(filepath, 'r') as f:  
        data = json.load(f)
    data_int = convert_index_to_int(data)
    for key, series in data_int.items():
        dfalldata["006F_DSSP_{}".format(key)] = pd.Series(series)
    log.info("Added .. {}".format(filename))

    #Functional Domains/Regions (UniProt)
    filename = '007A-StructuralProteome_UniprotFeatures_dict.json'
    filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
    with open(filepath, 'r') as f:  
        data = json.load(f)
    data_int = convert_index_to_int(data)
    for key, series in data_int.items():
        dfalldata["007A_FunctionalDomains_{}".format(key)] = pd.Series(series)
    log.info("Added .. {}".format(filename))

    #ALE mutants (index)
    filename =  '007B-StructuralProteome_ALE_mutations.json'
    filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
    with open(filepath, 'r') as f:  
        data = json.load(f)
    data_int = {}
    for k ,v in data.items():
        if str(k) in ['nan','NaN']:
            continue
        data_int.update({int(float(k)) : v})
    dfalldata['007B_LabMutant_ALE_csvIndex'] = pd.Series(data_int)
    log.info("Added .. {}".format(filename))

    #ALE mutants (grantham score)
    filename =   '007B-StructuralProteome_ALE_Grantham.json'
    filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
    with open(filepath, 'r') as f:  
        data = json.load(f)
    data_int = {}
    for k ,v in data.items():
        if str(k) in ['nan','NaN']:
            continue
        data_int.update({int(float(k)) : v})
    dfalldata['007B_LabMutant_ALE_VariantGS'] = pd.Series(data_int)
    log.info("Added .. {}".format(filename))

    #LTEE mutants (index)
    filename =  '007C-StructuralProteome_LTEE_mutations.json'
    filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
    with open(filepath, 'r') as f:  
        data = json.load(f)
    data_int = {}
    for k ,v in data.items():
        if str(k) in ['nan','NaN']:
            continue
        data_int.update({int(float(k)) : v})
    dfalldata['007C_LabMutant_LTEE_csvIndex'] = pd.Series(data_int)
    log.info("Added .. {}".format(filename))

    #LTEE mutants (grantham score)
    filename =   '007C-StructuralProteome_LTEE_Grantham.json'
    filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
    with open(filepath, 'r') as f:  
        data = json.load(f)
    data_int = {}
    for k ,v in data.items():
        if str(k) in ['nan','NaN']:
            continue
        data_int.update({int(float(k)) : v})
    dfalldata['007C_LabMutant_LTEE_VariantGS'] = pd.Series(data_int)
    log.info("Added .. {}".format(filename))

    
    
    #WT alleleome
    column_filenames = {"dfzIndex":"007D-StructuralProteome_global_dfz_index_WT_variation.json",
                        "Codon_DominantWTFreq":"007D-StructuralProteome_WT_CODON_dom_frequency.json",
                        "Codon_WTVariants":"007D-StructuralProteome_WT_CODON_seq_details.json",
                        "AA_DominantWTFreq":"007D-StructuralProteome_WT_AA_dom_frequency.json",
                        "AA_WTVariants":"007D-StructuralProteome_WT_AA_seq_details.json",
                       }
    for col, filename in column_filenames.items(): 
        filepath = op.join(qspaceDirs['DataOutput_dir'], filename)
        with open(filepath, 'r') as f:  
            data = json.load(f)
        data_int = {}
        for k ,v in data.items():
            if str(k) in ['nan','NaN']:
                continue
            data_int.update({int(float(k)) : v})
        dfalldata['007D_WT_alleleome_{}'.format(col)] = pd.Series(data_int)
        log.info("Added .. {}".format(filename))
    return dfalldata