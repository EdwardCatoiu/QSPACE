import os
import errno
import os.path as op
import json
from tqdm import tqdm_notebook as tqdm
import ast
import pandas as pd
import numpy as np
import urllib
# import urllib2
import requests
import copy

from Bio.PDB.Polypeptide import three_to_one

with open('qspace_directories.json','r') as f:
    qspaceDirs= json.load(f)
    
import logging
#import qspace directories
import sys
# sys.path.append('../')
# import qspace_directories as qspaceDirs

#import qspace functions
# sys.path.append(qspaceDirs.FunctionsDir)
import prepare_structures_for_BFS as prepBFS
import oligomerization

#import ssbio Functions
# sys.path.append(qspaceDirs.ssbioDir)
import ssbio
from ssbio.protein.sequence.seqprop import SeqProp
from ssbio.protein.structure.structprop import StructProp



updated_me_recipe = {'ABC-11-CPLX': {u'b0151': 2, u'b0152': 1, u'b0153': 1},
                     'ABC-13-CPLX':{u'b0654': 1, u'b0655': 1, u'b0652': 2, u'b0653': 1},
                     'ABC-18-CPLX': {u'b2148': 2, u'b2149': 1, u'b2150': 1},
                     'ABC-2-CPLX': {u'b1900': 1, u'b1901': 1, u'b4460': 2},
                     'ABC-26-CPLX': {u'b2677': 2, u'b2678': 2, u'b2679': 1},
                     'ABC-32-CPLX': {u'b0066': 2, u'b0067': 1, u'b0068': 1},
                     'ABC-33-CPLX': {u'b3566': 1, u'b3567': 1, u'b3568': 2},
                     'ABC-42-CPLX': {u'b4086': 2, u'b4087': 1, u'b4088': 1},
                     'ABC-46-CPLX': {u'b4227': 1, u'b4230': 1, u'b4231': 1, u'b4485': 1},
                     'ABC-49-CPLX': {u'b0829': 1, u'b0830': 1, u'b0831': 1, u'b0832': 1},
                     'ABC-6-CPLX': {u'b0886': 1, u'b0887': 1},
                     'ABC-53-CPLX': {u'b3201': 2, u'b4261': 1, u'b4262': 1},
                     'CPLX-165': {u'b1817': 1, u'b1818': 3, u'b1819': 3},
                     "EG11009-MONOMER"  : {u'b3035': 3},
                     "DIAMINOPIMDECARB-CPLX"  : {u'b2838': 2},
                       "ENTA-CPLX"  : {u'b0596': 4},
                     'DGTPTRIPHYDRO-CPLX' : {u'b0160': 6},
                     'CARBPSYN-CPLX' : {u'b0032': 4, u'b0033': 4},
                      'CPLX0-201':{u'b0825': 10},
                     u'EG11910-MONOMER_dimer_EG11911-MONOMER' :{u'b3951': 2,u'b3952': 1},
                     u'HYDROPEROXIDII-CPLX' :  {u'b1732': 4},
                     u'RNAPE-CPLX':{u'b3649': 1, u'b3295': 2, u'b3988': 1, u'b3987': 1, u'b2573': 1},
                     'Sec-CPLX': {'b0407': 1, 'b3175': 1, 'b3705': 1, 'b3300': 1, 'b0409': 1, 'b0408': 1, 'b3981': 1},
                     "NAP-CPLX":{u'b2203': 1, u'b2206': 1},
                     "NITRITREDUCT-CPLX" : {u'b3366': 1, u'b3365': 1},
                     "CPLX0-2081" : {u'b1200': 2, u'b1199': 2, u'b1198': 2},
                     "NRDACTMULTI-CPLX" :{u'b4238': 2,  u'b4237': 2},
                     "TMAOREDUCTI-CPLX" : {u'b0997': 1, u'b0996': 1},
                     "CPLX-153" : { u'b2415': 1},
                     "CPLX-168" : { u'b2417': 1,  u'b4240': 2},
                     "CPLX-157" :{ u'b2417': 1,  u'b1101': 2},
                     "CPLX0-231" : {u'b2092': 2, u'b2093': 1,  u'b2094': 1},
                     "EIISGA" : { u'b4193': 2, u'b4195': 1, u'b4194': 1},
                     "CPLX-154" :{u'b3722': 2},
                     "CPLX-166" :{u'b3599': 2},
                     "CPLX-167" : {u'b0679': 2},
                     "CPLX-163" : { u'b0731':  1},
                     "CPLX-158" : {u'b2167': 2, u'b2169': 1},
                      "LolCDE-CPLX" : {u'b1118': 1, u'b1116': 1, u'b1117': 2},
                     'CPLX0-3161':   {u'b2472': 2},
                     'CPLX0-7990': {'b3791': 2},
                     'CPLX0-263' :  {u'b0828': 2},
                     'APO-ENTB' :   {u'b0595': 2},
                     "CPLX0-7823":     {u'b1490': 2, u'b1489': 2},
                     "CPLX0-8218":     {u'b1490': 2},
                     "CPLX0-8199":     {u'b1489': 2},
                     "HMP-P-KIN-CPLX" : {u'b2103': 2},
                     "CPLX0-4" :  {u'b3240': 3, u'b3241': 6},
                     'METNIQ-METHIONINE-ABC-CPLX': {u'b0197': 1, u'b0199': 2, u'b0198': 2},
                     'ABC-35-CPLX' : {u'b2197': 2, u'b2198': 1, u'b2199': 1, u'b2200': 1, u'b2201': 2},
                     "MdtQ_trimer_efflux_pump" : {'b2139' : 3},
                     "MdtP_trimer_efflux_pump" : {'b4080' : 3},
                     "ABC-4-CPLX" : {u'b0864': 2, u'b0861': 1, u'b0860': 1, u'b0862': 1},
                     "ACYLCOASYN-CPLX":{u'b1805': 1},
                     "CPLX0-1924" : {u'b1252': 1, u'b3005': 2, u'b3006': 5, u'b3966': 1},
                     "CPLX0-1923_EG10155-MONOMER": {u'b1252': 1, u'b3006': 5, u'b2155': 1, u'b3005': 2},
                     'CPLX0-1923_EG10306-MONOMER' : {u'b1102': 1, u'b3006': 5, u'b3005': 2, u'b1252': 1},
                     'CPLX0-1941':{u'b1252': 1, u'b3006': 5, u'b0584': 1, u'b3005': 2},
                     'CPLX0-1942': {u'b0150': 1, u'b1252': 1, u'b3006': 5, u'b3005': 2},
                     'G6414-MONOMER_CPLX0-1923': {u'b0805': 1, u'b1252': 1, u'b3006': 5, u'b3005': 2},
                     'CPLX0-1943': {u'b4291': 1, u'b3006': 5, u'b3005': 2, u'b1252': 1},
                    
                    }