import os
import os.path as op
from tqdm import tqdm_notebook as tqdm
import pandas as pd
import json

import sys   
sys.path.append('functions')

from ssbio_functions import download_uniprot_file

with open('qspace_directories.json','r') as f:
    qspaceDirs= json.load(f)
    
import logging
