def find_real_chains(structure):
    res_list = ['ALA', 'ARG', 'ASN', 'ASP', 'GLN', 'CYS', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    chain_list = []
    for c in structure.structure.get_chains():
        found_real_chain = False
        for r in c.get_residues():
            if r.resname in ['DUM', 'HOH']:
                continue
            if r.resname in res_list and c.id != '_': #for swiss files 
                found_real_chain = True
                break
        if found_real_chain:
            chain_list.append(c.id)
    return chain_list