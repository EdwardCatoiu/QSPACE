import pandas as pd
import copy
from copy import deepcopy
from Bio import SeqIO
import os.path as op

import logging
log = logging.getLogger(__name__)


uniprot_total_feature_list = ['CARBOHYD','PEPTIDE','PROPEP','VARIANT','VAR_SEQ','INIT_MET','SIGNAL',
                             'PROPEP','TRANSIT','TOPO_DOM','TRANSMEM','INTRAMEM','DOMAIN','REPEAT','CA_BIND',
                             'ZN_FING','DNA_BIND','NP_BIND','REGION','COILED','MOTIF','COMPBIAS','ACT_SITE',
                             'METAL','BINDING','SITE','NON_STD','MOD_RES','LIPID','DISULFID','CROSSLNK','MUTAGEN',                             'UNSURE',
                             'CONFLICT','NON_CONS','NON_TER','HELIX','STRAND','TURN']

uniprot_secondary_struct_features = ['HELIX', 'TURN', 'STRAND', 'TOPO_DOM', 'TRANSMEM', "INTRAMEM"]

uniprot_mutagen_features = ['MUTAGEN','CONFLICT', 'VARIANT']

uniprot_active_site_features = ['CARBOHYD', 'CHAIN','PEPTIDE', 'PROPEP','VAR_SEQ', 'INIT_MET',
                                'SIGNAL', 'PROPEP','TRANSIT','DOMAIN','REPEAT', 'CA_BIND',
                                'ZN_FING', 'DNA_BIND','NP_BIND', 'REGION','COILED', 'MOTIF',
                                'COMPBIAS', 'ACT_SITE','METAL', 'BINDING','SITE', 'NON_STD', 
                                'MOD_RES', 'LIPID','DISULFID','CROSSLNK','UNSURE', 
                                'NON_CONS', 'NON_TER', 'INIT_MET']


def parse_uniprot_features_from_txt_jan2020(uniprot_file,selected_features= uniprot_total_feature_list):
    
    
    """ Read the uniprot file and parse it for features (FT)

        Args:
            uniprot_file : <str> path to uniprot file
            
        Returns:
            feature_sites: <pd.Dataframe> dataframe of uniprot features, seq start/end, descriptions
            features_available: <list> all unique features from uniprot
            """
        
        
        

    feature_sites = pd.DataFrame(columns = ['feature', 'seq_start', 'seq_end', 'description'])
    index = 0
    with open(uniprot_file, 'r') as txt:

        cache_txt = txt.read()

        for l in cache_txt.splitlines():
            use_previous_line_info = False
            feature_tag = l[:5].strip()

            line = l[5:].strip()
            words = line.split()
            
            if feature_tag != 'FT':
                continue
            if words[0] =='CHAIN':
                continue
            if words[0] not in selected_features :
                use_previous_line_info = True
            #print l, words, use_previous_line_info
            if not use_previous_line_info:
                if '?' in words[1:]:
                    continue
                
                if len(words) == 3:
                    [ft, start,end] = words
                else:
                    ft = words[0]
                    words[1].split('.')
                    
                    start = words[1].split('.')[0]
                    end = words[1].split('.')[-1]
                try:
                    feature_sites.loc[index, 'feature'] = ft
                    feature_sites.loc[index, 'seq_start'] = int(start)
                    feature_sites.loc[index, 'seq_end'] = int(end)
                except ValueError:
                    log.info("Could not find Feature Start/End Residue | {} {} {} {} ".format(op.basename(uniprot_file),ft, start, end))
                description = ''
                for w in words[3:]:
                    if 'ECO' in w and 'Rule' not in w:
                        continue
                        
                    elif 'ECO' in w:
                        w = w.split('|')
                        for e in w:
                            if 'Rule' in e:
                                description = description + ' ' + e
                    else:
                        description = description + ' ' + w

                feature_sites.loc[index , 'description'] = description
                index = index + 1
            else:
                
                index = index - 1
                if index == -1: 
                    index = 0 
                    continue
                description = feature_sites.loc[index , 'description']
                for w in words:
                    if 'ECO' in w and 'Rule' not in w:
                        continue
                        
                    elif 'ECO' in w:
                        w = w.split('|')
                        for e in w:
                            if 'Rule' in e:
                                description = description + ' ' + e
                    else:
                        description = description + ' ' + w


                feature_sites.loc[index, 'description'] = description
                index = index + 1


    features_available = []
    for index, row in feature_sites.iterrows():


        if row.feature =='CHAIN':
            continue

        if row.feature in ['STRAND', 'HELIX', 'TURN' ] or row.description.strip('') == '':
            feature = row.feature

        else:
            feature = row.feature + '_info_' + row.description

        if feature not in features_available:
            features_available.append(feature)

    #features_available =  feature_sites.feature.drop_duplicates().tolist()
    if 'CHAIN' in features_available:
        features_available.remove('CHAIN')
    return feature_sites, features_available

def parse_uniprot_features_from_txt(uniprot_file):
    
    
    """ Read the uniprot file and parse it for features (FT)

        Args:
            uniprot_file : <str> path to uniprot file
            
        Returns:
            feature_sites: <pd.Dataframe> dataframe of uniprot features, seq start/end, descriptions
            features_available: <list> all unique features from uniprot
            """
        
        
        

    feature_sites = pd.DataFrame(columns = ['feature', 'seq_start', 'seq_end', 'description'])
    index = 0
    with open(uniprot_file, 'r') as txt:

        cache_txt = txt.read()

        for l in cache_txt.splitlines():
            use_previous_line_info = False
            feature_tag = l[:5].strip()

            line = l[5:].strip()
            words = line.split()
            
            if feature_tag != 'FT':
                continue
            if words[0] not in uniprot_total_feature_list:
                use_previous_line_info = True

            if not use_previous_line_info:
                if '?' in [words[1], words[2]]:
                    continue
                
                
                feature_sites.loc[index, 'feature'] = words[0]
                feature_sites.loc[index, 'seq_start'] = int(words[1])
                feature_sites.loc[index, 'seq_end'] = int(words[2])

                description = ''
                for w in words[3:]:
                    if 'ECO' in w and 'Rule' not in w:
                        continue
                        
                    elif 'ECO' in w:
                        w = w.split('|')
                        for e in w:
                            if 'Rule' in e:
                                description = description + ' ' + e
                    else:
                        description = description + ' ' + w

                feature_sites.loc[index , 'description'] = description
                index = index + 1
            else:
                index = index - 1
                description = feature_sites.loc[index , 'description']
                for w in words:
                    if 'ECO' in w and 'Rule' not in w:
                        continue
                        
                    elif 'ECO' in w:
                        w = w.split('|')
                        for e in w:
                            if 'Rule' in e:
                                description = description + ' ' + e
                    else:
                        description = description + ' ' + w


                feature_sites.loc[index, 'description'] = description
                index = index + 1


    features_available = []
    for index, row in feature_sites.iterrows():


        if row.feature =='CHAIN':
            continue

        if row.feature in ['STRAND', 'HELIX', 'TURN' ] or row.description.strip('') == '':
            feature = row.feature

        else:
            feature = row.feature + '_info_' + row.description

        if feature not in features_available:
            features_available.append(feature)

    #features_available =  feature_sites.feature.drop_duplicates().tolist()
    if 'CHAIN' in features_available:
        features_available.remove('CHAIN')
    return feature_sites, features_available


def extract_features_from_uniprot_df(df, feature_list):
    """ Convert the uniprot dataframe into a dictionary

    Args:
        df : <pd.Dataframe> uniprot feature dataframe
        feature_list: <list> all unique features from uniprot

    Returns:
        uniprot_feature_dict: <dict> {uniprot_feature : [ [ res1, res2, res3...] ,  [ res15, res16, res17...] ...]
            the values are seperated by their occurence in the uniprot metadata
    """
    
    
    
    uniprot_feature_dict = {}
    for feature in feature_list:
        addtl_info = False
        if '_info_' in feature:
            addtl_info = True

        seq_list = []
        for index, row in df.iterrows():
            if not addtl_info:
                if row.description != '' or row.feature != feature:
                    continue
            elif row.feature not in feature:
                continue
            elif row.description == '' or row.description not in feature:
                continue
#             print index, feature, '\t', '---->',row.feature, row.description, str(addtl_info)
            start = row.seq_start
            end = row.seq_end
            if str(start) =='nan' or str(end) =='nan':
                log.info("Could not map Feature Start/End Residue | {} {} {} {} ".format(index,feature, start, end))
                continue
            if 'CROSSLNK' in feature or 'DISULFID' in feature:
                ft_range = [start, end]
            else:
                ft_range = range(start, end+1)
            seq_list.append(list(set(ft_range)))
        uniprot_feature_dict.update({feature : seq_list})
    return uniprot_feature_dict


def combine_crosslink_if_same(feature_list, feature_dict):
    """ if "CROSSLNK" is a feature in the uniprot metadata, then we need to parse it a different way to ensure only the seq
        start and ends make it to the uniprot feature dictionary.

    Args:
        feature_dict : <dict> {uniprot_feature : [ [ res1, res2, res3...] ,  [ res15, res16, res17...] ...]
        feature_list: <list> all unique features from uniprot

    Returns:
        new_feature_dict: <dict> {uniprot_feature : [ [ res1, res2, res3...] ,  [ res15, res16, res17...] ...]
            the values are seperated by their occurence in the uniprot metadata, CROSSLINK features are combined
            if they are the same
    """
    
    
    #this code combines cross link features (such a pain in the ass)
    
    new_feature_dict = copy.deepcopy(feature_dict)
    
    crosslink_list = []
    crosslink_dict = {}
    for feature in feature_list:
        if 'CROSSLNK' in feature:
#             print feature, feature_dict[feature]
            crosslink_list.append(feature)
            crosslink_dict.update({feature : feature_dict[feature]})
            del new_feature_dict[feature]
    
    residue_list = set()
    common_residues = set()
    for crslink, occurence_list in crosslink_dict.items():
        for single_crosslink_occurence, single_crosslink_residue_list  in enumerate(occurence_list):
#             print single_crosslink_occurence, single_crosslink_residue_list
#             print residue_list, '<<<'
#             print residue_list.intersection(set(single_crosslink_residue_list))
            if residue_list.intersection(set(single_crosslink_residue_list)) == set():
                residue_list = residue_list.union(set(single_crosslink_residue_list))
                
            else: 
                common_residues = common_residues.union(residue_list.intersection(set(single_crosslink_residue_list)))
                
    
    to_use_crosslinks = copy.deepcopy(crosslink_dict)
    
    shared_residue_dict = {}
    for shared_residue in common_residues:
        for crslink, crslink_residue_lists in crosslink_dict.items():
            for residue_list in crslink_residue_lists:
                if shared_residue in residue_list:
                    if shared_residue in shared_residue_dict.keys():
                        shared_residue_dict[shared_residue] += [crslink]
                    else:
                        shared_residue_dict.update({shared_residue: [crslink]})
                    
#                     print shared_residue, crslink
    
    new_crosslink_dict = {}
    for shared_residue, identical_crosslinks in shared_residue_dict.items():
        new_crosslink_list = []
        new_crosslink_name = ''
        
        for crslink in identical_crosslinks:
            new_crosslink_name = new_crosslink_name +  '(%s)' %crslink.split('(')[1].split(')')[0]
            
            for i, residues in enumerate(crosslink_dict[crslink]):
                for r in residues:
                    if r not in new_crosslink_list:
                        new_crosslink_list.append(r)
                        
                        
        new_crosslink_name = 'CROSSLNK_info_ ' + new_crosslink_name                
#         print '>>>>' , shared_residue , new_crosslink_list, new_crosslink_name
        new_crosslink_dict.update({new_crosslink_name : [new_crosslink_list]})
    
    new_feature_dict.update(new_crosslink_dict)
    
    return new_feature_dict


def make_nextGenome_feature_dataframe(gene_name, uniprot_text_path, uniprot_fasta_path):
    
    
    if uniprot_text_path is not None:
        
        uniprot_df, feature_list = parse_uniprot_features_from_txt_jan2020(uniprot_text_path)
        feature_dict = extract_features_from_uniprot_df(df=uniprot_df, feature_list= feature_list)
        feature_dict = combine_crosslink_if_same(feature_list, feature_dict)
    else:
        feature_dict = {}
    
    res_dict = {}
    aa_seq = SeqIO.read(open(uniprot_fasta_path), 'fasta')
    for i , AA in enumerate(str(aa_seq.seq)):
        res_dict.update({"{}_{}".format(gene_name,i + 1 )  : {'AA_Seq' : AA, 
                                                              'AA_SeqNum' : i + 1, 
                                                              'Gene_ID' : gene_name,
                                                              'UniProt' : op.basename(uniprot_fasta_path) }})

    
    
    
    for f, reslist_all in feature_dict.items():
        if 'info' in f:
            [domain, metadata] = f.split('_info_ ')
        else:
            domain = f
            if domain in ['HELIX', 'STRAND', 'TURN']:
                metadata = domain
            else:
                metadata = 'Unknown'
            
        for reslist in reslist_all:
            for resnum in reslist:
                index = "{}_{}".format(gene_name, resnum)
                if index in res_dict.keys():
                    if domain in res_dict[index]:
                        if metadata in res_dict[index][domain]:
                            continue
                        else:
                            res_dict[index][domain] += [metadata]
                    else:
                        res_dict[index].update({domain : [metadata]})
                else:
                    
                    res_dict.update({index : {domain : [metadata]} })
                 
        
                        
            
    testdf = pd.DataFrame.from_dict(res_dict)
    return testdf

