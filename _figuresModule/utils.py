import os
import os.path as op
import json
from tqdm import tqdm_notebook as tqdm
import ast
import pandas as pd
import numpy as np
import copy

from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

with open('qspace_directories.json','r') as f:
    qspaceDirs= json.load(f)
    
import sys
# sys.path.append('../')
# import qspace_directories as qspaceDirs

#import qspace functions
sys.path.append(qspaceDirs['FunctionsDir'])
import venn
from venn import *

import seaborn as sns
