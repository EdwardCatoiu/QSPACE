import os 
import os.path as op
import pandas as pd
import numpy as np
import seaborn as sns

###ssbio imports
import sys
from ssbio.protein.structure.utils.structureio import StructureIO

### qspace imports
import protein_geometry



AA_resnames = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}.keys()


def find_coords(s):

    coorindate_list = {}
    for c in s.get_chains():
        for r in c.get_residues():
            if r.resname not in AA_resnames:
                continue
            
            try:
    #             print r.child_dict['CA'].coord
                coorindate_list.update({c.id + '_' + str(r.id[1]) : list(r.child_dict['CA'].coord)})
            except KeyError:
                continue
#     print len(coorindate_list)
    return coorindate_list

def color_string_structure(df_atoms, inner_fence, color_count = 0, view = False):
  
    
    #NEWWWWWW
    current_palette = sns.color_palette("hls", 8).as_hex()  
    if color_count >= len(current_palette):
        color_count = 0
    # df_atoms = df_atoms.sort_values(by = 'distance', ascending=False)
    remove_index = []
    for index, row in df_atoms.iterrows():
        if row.distance > inner_fence:
            remove_index.append(index)
            if not view:
                continue
            [chain,r] = index.split('_')
            view.add_ball_and_stick(selection = ':{} and (not hydrogen and {})'.format(chain, r), 
                                        scale = 5, opacity = 0.25, color = current_palette[color_count])
        if row.distance == min(df_atoms.distance):
#             print index
            if not view:
                continue
            [chain,r] = index.split('_')
            view.add_ball_and_stick(selection = ':{} and (not hydrogen and {})'.format(chain, r), 
                                        scale = 15, opacity = 0.6, color = current_palette[color_count])
            
    df_atoms= df_atoms[df_atoms.index.isin(remove_index) == False]
    p0 = [np.mean(df_atoms[0]), np.mean(df_atoms[1]),np.mean(df_atoms[2])]

    color_count += 1
    
    return df_atoms, inner_fence,view,remove_index, color_count

def visualize_outliers(filepath, view = False):
    
    if not op.exists(filepath):
        return IOError
    
    def get_df_atoms_fence(df_atoms):
        p0 = [np.mean(df_atoms[0]), np.mean(df_atoms[1]),np.mean(df_atoms[2])]
        for index in df_atoms.index:
            d = protein_geometry.distance_two_points(p0, coords[index])
            df_atoms.loc[index, 'distance'] = d
        q75, q25 = np.percentile(df_atoms['distance'], [75 ,25])
        iqr = q75 - q25
        inner_fence = 3.0 * iqr + q75
        return df_atoms, inner_fence
    
    
     
    if view:
        view = nglview.show_structure_file(filepath)
        view.clear_representations()
        view.add_cartoon(selection='protein', color = 'silver')
    s = StructureIO(filepath)
    p = s.structure[0]
    coords = find_coords(p)
    df_atoms = pd.DataFrame.from_dict(coords, orient = 'index')
#     print len(df_atoms), 'len'
#     return df_atoms, coords
    color_count = 0
    remove_index = ['start']
    total_residues_removed = []
    inner_fence= 0 
    total_residues = len(coords)
    iteration_count = 0 
#     print len(df_atoms), 
    if len(df_atoms) > 0:
        while len(remove_index) > 0:
            df_atoms, inner_fence = get_df_atoms_fence(df_atoms)
            df_atoms, inner_fence,view,remove_index , color_count= color_string_structure(df_atoms = df_atoms, inner_fence= inner_fence , view= view,  color_count = color_count)
            total_residues_removed += remove_index
            iteration_count += 1
    if len(total_residues_removed)> 0:
        pass
#         print len(total_residues_removed), '\t', op.basename(filepath)
    return df_atoms,inner_fence, view,total_residues_removed,total_residues ,iteration_count



def atom_line(atom):
    line = '{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2}\n'.format('ATOM',
                                                                                                                       atom.serial_number,
                                                                                                                       atom.name,
                                                                                                                       '',
                                                                                                                       atom.parent.resname,
                                                                                                                       'X',
                                                                                                                       atom.parent.id[1],
                                                                                                                       '',
                                                                                                                       atom.coord[0],
                                                                                                                       atom.coord[1],
                                                                                                                       atom.coord[2],
                                                                                                                       atom.occupancy,
                                                                                                                       atom.bfactor,
                                                                                                                       atom.element)
    
    return line
    