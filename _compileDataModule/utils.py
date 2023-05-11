import pandas as pd
import os
import os.path as op
from tqdm import tqdm_notebook as tqdm
import json
import copy

import sys
with open('qspace_directories.json','r') as f:
    qspaceDirs= json.load(f)
    
sys.path.append(qspaceDirs['FunctionsDir'])
import logging
