import pandas as pd
from tqdm import tqdm_notebook as tqdm
import os
import os.path as op
import json
import shutil
import sys
#import qspace directories
import sys
# sys.path.append('../')
# import qspace_directories as qspaceDirs

#import qspace functions

with open('qspace_directories.json','r') as f:
    qspaceDirs= json.load(f)
sys.path.append(qspaceDirs['FunctionsDir'])
import read_uniprot_text
import uniprot_features
# import prepare_structures_for_BFS as prepBFS
# import oligomerization
# import bfs_utils

#import ssbio Functions
from ssbio.protein.structure.structprop import StructProp

###### webserver OPM
from selenium.webdriver.chrome.options import Options
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
import time

import requests
from bs4 import BeautifulSoup


def scrape_orient(table_orientation):
    data = table_orientation.find_all('tr', attrs={"class":["row light"]})
    datalist = []
    for value in list(data[0]):
        try:
            value = value.contents
        except AttributeError:
            continue
#         print "---{}---".format(value)
        value = value[0].replace(' ','')
        value = value.split(u'\xb1')[0]
        value = value.split(u'kcal')[0]
        value = float(value)
        datalist += [value]
    return {"Depth/Hydrophobic Thickness" : datalist[0],
            "delta_G" : datalist[1],
            "Tilt angle" : datalist[2],
           }

def scrape_residues(table_embedded):
    tbody = table_embedded.find_all('tbody')[0]

    embedded_residues= False
    transmembrane_segments= False
    EmbeddedResiduesDict = {}
    TransmembraneSegmentsDict = {}
    for row in tbody.find_all('tr'):
        if row.attrs == {'class': ['row', 'dark']}:
            if 'Embedded_residues' in row.getText():
                embedded_residues = True
                transmembrane_segments = False

                continue
            if 'Transmembrane_secondary_structure_segments' in row.getText():
                embedded_residues = False
                transmembrane_segments = True
                continue

        if row.attrs == {'class': ['row', 'light']}:
            if embedded_residues:
                [chaininfo,tilt,residues] = row.find_all('td')
    #                 chain = chaininfo[0].getText()
    #             print chaininfo.getText(), tilt.getText(),residues.getText().split(',')
                chainText = chaininfo.getText()
                chainText = chainText.replace(' ','')
                try:
                    tiltText = float(tilt.getText())
                except ValueError:
                    tiltText = tilt.getText()

                residuesList = residues.getText().split(',')
                for segment in residuesList:
                    try:
                        if chainText in EmbeddedResiduesDict:
                            EmbeddedResiduesDict[chainText] += [int(segment)]
                        else:
                            EmbeddedResiduesDict.update({chainText : [int(segment)]})

                    except ValueError:
                        segment = segment.split('-')
                        segment_list = list(range(int(segment[0]),int(segment[-1])+1))
                        if chainText in EmbeddedResiduesDict:
                            EmbeddedResiduesDict[chainText] += segment_list
                        else:
                            EmbeddedResiduesDict.update({chainText : segment_list})

            elif transmembrane_segments:
                [chaininfo,tilt,residues] = row.find_all('td')
    #                 chain = chaininfo[0].getText()
    #             print chaininfo.getText(), tilt.getText(),residues.getText().split(',')
                chainText = chaininfo.getText()
                chainText = chainText.replace(' ','')
                try:
                    tiltText = float(tilt.getText())
                except ValueError:
                    tiltText = tilt.getText()
                residuesList = residues.getText().split(',')
                residuesInSegment = []
                for segment in residuesList:
                    segment = segment.split('(')[-1]   
                    segmentStart = int(segment.split('-')[0])
                    segmentEnd = int(segment.split(')')[0].split('-')[-1])
    #                 print segment, segmentStart,segmentEnd
                    residuesInSegment.append(list(range(segmentStart,segmentEnd+1)))

                if chainText in TransmembraneSegmentsDict:
                    TransmembraneSegmentsDict[chainText] += residuesInSegment
                else:
                    TransmembraneSegmentsDict.update({chainText : residuesInSegment})
    return EmbeddedResiduesDict,    TransmembraneSegmentsDict

def getPDBlink(table_output):

    t = table_output.find_all('tbody')[0]
    t = t.find_all('td')[0]
    pdbLink = str(t).split('href="')[1].split('.pdb')[0] + '.pdb'
    return pdbLink

def opmNeeded(query_structureIds,
              OPMstructureFolder = qspaceDirs['opmOutputStructuresDir'],
              force_rerun=False):
    
    NotSolvedInfile = op.join(qspaceDirs['DataOutput_dir'], '005B-OPM-not_solved_dict.json')
    if op.exists(NotSolvedInfile):
        with open(NotSolvedInfile, 'rb') as f:
            not_solved_dict = json.load(f)
        
    
    opm_needed = []
    for s_id in query_structureIds:
        if op.exists(op.join(OPMstructureFolder, "OPMstructure-{}out.pdb".format(s_id.replace('-','_')))) and not force_rerun:
            if s_id in not_solved_dict:
                del not_solved_dict[s_id]
                with open(NotSolvedInfile, 'wb') as f:
                    json.dump(not_solved_dict,f)
            
            continue
        opm_needed += [s_id]
        
    print ("{} of {} still need OPM".format(len(opm_needed), len(query_structureIds)))
    return opm_needed
