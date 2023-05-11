import pandas as pd
import os
import os.path as op
from tqdm import tqdm_notebook as tqdm
import json
import copy

import sys
#import qspace directories
with open('qspace_directories.json','r') as f:
    qspaceDirs= json.load(f)
    
sys.path.append(qspaceDirs['FunctionsDir'])
import uniprot_features
import logging


#import ssbio Functions
from ssbio.protein.sequence.properties.residues import  grantham_score
import ssbio.protein.sequence.utils.alignment as ssbioaln