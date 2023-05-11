import os
import os.path as op
import json
import pandas as pd
import numpy as np


def find_GO(filepath):
    go_list_process = {}
    go_list_function = {}
    go_list_component = {}
    with open(filepath, 'r') as txt:

        cache_txt = txt.read()

        for l in cache_txt.splitlines():
            use_previous_line_info = False
            feature_tag = l[:5].strip()

            line = l[5:].strip()
            words = line.split()

            if feature_tag != 'DR':
                continue
            if words[0] !='GO;':
                continue
#             print words
            if words[2][0] == 'P':
                go_list_process.update({words[1].rstrip(';') : ' '.join(words[2:-1]).split(':')[1].rstrip(';') })
                
            if words[2][0] == 'C':
                go_list_component.update({words[1].rstrip(';'): ' '.join(words[2:-1]).split(':')[1].rstrip(';')})
            if words[2][0] == 'F':
                go_list_function.update({words[1].rstrip(';'): ' '.join(words[2:-1]).split(':')[1].rstrip(';')})
    
    txt.close()
    return go_list_process,go_list_function, go_list_component

