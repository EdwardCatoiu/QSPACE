import pandas as pd
import copy

def parse_uniprot_TM_features_from_txt(uniprot_file):
    """
    args : str (infile path to uniprot file)
    parses the Uniprot and returns a dataframe with just "FT" rows
    eddies version
    """
    
    uniprot_membrane_sites = pd.DataFrame(columns = ['type', 'seq_start', 'seq_end',
                                                         'location','location_transmem_start',
                                                         'location_transmem_end'])
    index = 0
    with open(uniprot_file, 'r') as txt:
       
        cache_txt = txt.read()
     
        for l in cache_txt.splitlines():
#             print l
            test_tag = l[:5].strip()
            line = l[5:].strip()
            line = line.replace('.',' ')
            words = line.split()
#             print '\t',words
            if test_tag != 'FT':    
#                 if '/note' not in str(words):
                continue
                    
#             print '\t-->{}<--'.format(words)
    
            if words[0] not in ['TOPO_DOM', 'TRANSMEM', 'INTRAMEM']:
                if "/note=" in words[0] and index >0:
                    location = words[0].split('=')[1]
                    location = location.replace('"', '')
                    location = location.replace("'", '')
                    if location not in ['Cytoplasmic','Periplasmic','Extracellular','Helical']:
                        continue
                    uniprot_membrane_sites.loc[index-1, 'location'] = location
#                     print l
                continue
            uniprot_membrane_sites.loc[index, 'type'] = words[0]
            uniprot_membrane_sites.loc[index, 'seq_start'] = words[1]
            try:
                uniprot_membrane_sites.loc[index, 'seq_end'] = words[2]
            except IndexError:
                uniprot_membrane_sites.loc[index, 'seq_end'] = words[1]
            if len(words) >= 4:
                uniprot_membrane_sites.loc[index, 'location'] = words[3].split(';')[0].split('.')[0]
            index = index + 1

            
    df = copy.deepcopy(uniprot_membrane_sites)       
    for index, row in df.iterrows():
        if index == max(df.index) or index == min(df.index):
            #print 'FIRST/LAST ENTRY', row.seq_start, row.seq_end
            continue

        #identify previous and next features
        prev_index = index - 1
        next_index = index + 1
        prev_seq_end = df.loc[prev_index, 'seq_end']
        next_seq_start = df.loc[next_index, 'seq_start']
        prev_location = df.loc[prev_index, 'location']
        next_location = df.loc[next_index, 'location']
        if prev_seq_end == '?' or next_seq_start == '?':
            continue

        #previous domain is topological
        if row.type in  ['TRANSMEM', 'INTRAMEM'] and df.loc[prev_index, 'type'] == 'TOPO_DOM':
            #print index
            #print  prev_location,'\t', prev_seq_end, '>>>>>>>>>',row.seq_start,'\t', row.type
            df.loc[index, 'location_transmem_start'] = prev_location

        #next domain is topological
        if row.type in  ['TRANSMEM', 'INTRAMEM'] and df.loc[next_index, 'type'] == 'TOPO_DOM':
            #print row.type,'\t', row.seq_end, '>>>>>>>>>', next_seq_start, '\t',next_location
            df.loc[index, 'location_transmem_end'] = next_location
            
            
    def get_features_of_membrane_annotation(df):    #self.TM_features
        feature_list = []
        #tm and topo query
        just_topo_residues = df[df['type'].str.contains('TOPO_DOM')]
        just_tm_residues = df[df['type'].str.contains('MEM')]
        ###############################################
        #feature 3: no tm domains from uniprot
        if len(just_tm_residues) == 0:
            feature_list.append('No_TM_Domains')
            return feature_list
        ###############################################
        #feature 4: no topological information, some TM domains
        if len(just_topo_residues) == 0:
            feature_list.append('No_Topo_Domains')
            return feature_list
        ###############################################
        first_index = min(df.index)
        last_index = max(df.index)
        #feature 5 : TM is first!
        if 'TOPO_DOM' not in df.loc[first_index, 'type']:
            feature_list.append('TM_domain_first')
        #feature 6 : TM is last!
        if 'TOPO_DOM' not in df.loc[last_index, 'type']:
            feature_list.append('TM_domain_last')
        ###############################################    
        #feature 7: multiple membrane domains in a row!
        increment = 1
        topo_index_list = []
        for index_m_u, row_m_u in df.iterrows():

            if row_m_u.type == 'TOPO_DOM':
                #print index_m_u
                topo_index_list.append(index_m_u)

        count = 0 
        feature_count = 0 
        for topo_index in topo_index_list:
            if topo_index == topo_index_list[0]:
                continue
            else:
                index_difference = topo_index - topo_index_list[count]
            if index_difference > 2:
                feature_count = feature_count + 1
            count = count + 1
        if feature_count > 0:
            feature_list.append('Multiple_BackToBack_TMs_{}_times'.format(feature_count))
        return feature_list
    
    TM_features = get_features_of_membrane_annotation(df)
    #######################
    
    #######################
    #  create_uniprot_tm_dataframe
    if 'No_Topo_Domains' in TM_features:
        domain_residues_1 = []
        domain_residues_2 = []
        count = 0
        temp_dict ={}
        for index, row in df.iterrows():
            #find odd / even count, regardless of index
            if count % 2 == 0 :
                domain_residues_1.append(row.seq_start)
                df.loc[index, 'location_transmem_start'] = 'TM_1'
                domain_residues_2.append(row.seq_end)
                df.loc[index, 'location_transmem_end'] = 'TM_2'

            elif count % 2 == 1:
                domain_residues_2.append(row.seq_start)
                df.loc[index, 'location_transmem_start'] = 'TM_2'
                domain_residues_1.append(row.seq_end)
                df.loc[index, 'location_transmem_end'] = 'TM_1'
            count = count + 1
        temp_dict.update({'TM_1' : domain_residues_1})
        temp_dict.update({'TM_2' : domain_residues_2})

    elif 'TM_domain_first' in TM_features or 'TM_domain_last' in TM_features:


        all_locations = set(df.location_transmem_start.unique().tolist() + df.location_transmem_end.unique().tolist())
        leaf_locations = copy.deepcopy(all_locations)
        for loc in all_locations:
            if loc is None or str(loc) == 'nan':
                leaf_locations.remove(loc)

        use_remaining_loc = False
        if len(leaf_locations) == 1 :
            only_loc = leaf_locations.pop()
            if only_loc =='Cytoplasmic':
                use_remaining_loc = True
                remaining_loc = 'Periplasmic_C'
            if only_loc =='Extracelluar':
                use_remaining_loc = True
                remaining_loc = 'Periplasmic_E'

        #assigns missing locations to the remaining location (periplasm) for when only Extracellular or Cytoplasmic locations 
        #are present. only works if just_tm_residues is VERY SHORT (1 row)

        first_index = min(df.index)
        second_index = first_index + 1

        last_index  = max(df.index)
        semi_last_index = last_index - 1

        index_dict = {first_index : [second_index, 'location_transmem_end', 'location_transmem_start' ],
                     last_index : [semi_last_index, 'location_transmem_start','location_transmem_end'] }



        for feature_index, adjacent_index in index_dict.items():
            #print feature_index

            #searches for TM residues in first/last position
            if df.loc[feature_index, 'type'] not in ['TRANSMEM', 'INTRAMEM']:
                continue

            if df.loc[adjacent_index[0], 'type'] == 'TOPO_DOM':
                topo_location = df.loc[adjacent_index[0], 'location']
                df.loc[feature_index, adjacent_index[1]] = topo_location

            #skips if the distance between start/end is less than 15 (not certain its a true TM)
            if int(df.loc[feature_index,'seq_end']) - int(df.loc[feature_index,'seq_start']) < 15:
                continue

            #if there is only one leaflet in the uniprot and its CYTO OR EXTRA
            if use_remaining_loc == True:
                df.loc[feature_index, adjacent_index[1]] = remaining_loc
                print ('remaining_loc activated')

            else:
                #print leaf_locations
                for loc in leaf_locations:
                    if df.loc[feature_index, adjacent_index[1]] != loc:
                        df.loc[feature_index, adjacent_index[2]] = loc
    return df, TM_features


def extract_residues_from_uniprot_TMdataframe(df,full_helix = False, proximal_distance = 3):
    if full_helix == False:
        u_leaf_dict = {}
        all_tm_locations = set(df.location_transmem_end.drop_duplicates().tolist() + df.location_transmem_start.drop_duplicates().tolist())
        all_tm_locations.update(['Unknown_Residues'])#finds the unique locations + The unknown ones
        for loc in all_tm_locations:
            if type(loc) != str:
                continue
            domain_residues = []
            for index, row in df.iterrows():
                if row.type == 'TOPO_DOM':
                    continue
                #now we assign TM/IM residues to each location
                #we update a dict as well
                if loc != 'Unknown_Residues':
                    if row.location_transmem_start == loc:
                        domain_residues.append(int(row.seq_start))
                    if row.location_transmem_end == loc:
                        domain_residues.append(int(row.seq_end))
                else:
                    if str(row.location_transmem_start) in [None, 'nan', 'None'] and row.seq_start != '?':
                        domain_residues.append(int(row.seq_start))
                    if str(row.location_transmem_end) in [None, 'nan', 'None'] and row.seq_end != '?':
                        domain_residues.append(int(row.seq_end))
            u_leaf_dict.update({loc : domain_residues})


        if 'Extracellular' in u_leaf_dict.keys() and 'Periplasmic' in u_leaf_dict.keys():
            if 'Periplasmic_E' not in u_leaf_dict.keys():
                u_leaf_dict.update({'Periplasmic_E' : u_leaf_dict['Periplasmic']})
            else:
                u_leaf_dict['Periplasmic_E'] += u_leaf_dict['Periplasmic']
            del u_leaf_dict['Periplasmic']

        if 'Cytoplasmic' in u_leaf_dict.keys() and 'Periplasmic' in u_leaf_dict.keys():
            if 'Periplasmic_C' not in u_leaf_dict.keys():
                u_leaf_dict.update({'Periplasmic_C' : u_leaf_dict['Periplasmic']})
            else:
                u_leaf_dict['Periplasmic_C'] += u_leaf_dict['Periplasmic']
            del u_leaf_dict['Periplasmic']

        return u_leaf_dict

    elif full_helix == True:
        TM_dict = {}
        count = 0 
        for index, row in df.iterrows():
            if row.type != 'TRANSMEM':
                continue
            if row.seq_start == '?' or row.seq_end == '?':
                continue
            count = count + 1
            transmem_key = row.type + '_%i' %count    

            full_tm_list = range(int(row.seq_start) + proximal_distance, int(row.seq_end) + 1 - proximal_distance)

            TM_dict.update({transmem_key : full_tm_list})

        return TM_dict