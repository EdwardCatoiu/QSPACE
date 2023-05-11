import pandas as pd
from tqdm import tqdm_notebook as tqdm

def get_PDB_gene_maps(mmSeqs,dfpdb_mapped):

    pdbMaps = {}
    for blattner, blattner_maps in tqdm(mmSeqs.items()):
        for pdb_entry, entity_maps in blattner_maps.items():
            dfp = dfpdb_mapped[dfpdb_mapped.pdb_entry == pdb_entry]

            for entity_id, scoreDict in entity_maps.items():
    #             dfp_es = dfp[dfp.entity_id == int(entity_id)]
    #             asym_ids = ast.literal_eval(dfp_es.asym_ids.tolist()[0])
                score = scoreDict['score']
        
                if pdb_entry not in pdbMaps:
                    pdbMaps.update({pdb_entry : {blattner :  {entity_id : score} }})
                elif blattner not in pdbMaps[pdb_entry]:
                    pdbMaps[pdb_entry].update({blattner : {entity_id :  score} })
                elif entity_id not in pdbMaps[pdb_entry][blattner]:
                    pdbMaps[pdb_entry][blattner].update({entity_id : score})
    return pdbMaps