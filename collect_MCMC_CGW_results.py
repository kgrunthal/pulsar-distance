# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 13:41:10 2024

@author: kgrun
"""

import json
import glob, os

import numpy as np


basepath='/u/kgrunthal/HD/'
n_params = 8

for lmc in [8.5, 9.0, 9.5]:
    for pd in [1.0, 1.5, 2.0]:
        paths = sorted(glob(basepath+'/MCMCout_CGW{}_pd{}_*'.format(lmc, pd)))
        iterations = len(paths)
        
        # create outputfile
        f = open(basepath+'/parameters_lmc{}_pd{}.txt'.format(lmc, pd), 'a')
        
        #iterate over realisations
        for i, path in enumerate(paths):
            
            with open(path, 'r') as current_file:
                maxlike_dict = json.load(current_file)
            current_file.close()
            
            outline = np.zeros(len(3*n_params))
            
            if i == 0:
                headline =''
                for j, key in enumerate(maxlike_dict.keys()):
                    if j == 0:
                        headline += '{}\t{}_ll\t{}_ul'.format(key, key, key)
                    else:
                        headline += '\t{}\t{}_ll\t{}_ul'.format(key, key, key)
                    
                    pos = int(3*j)
                    outline[pos] = maxlike_dict[key][1]
                    outline[pos+1] = maxlike_dict[key][0]
                    outline[pos+2] = maxlike_dict[key][2]
                    
                    f.write(headline + '\n')
                    f.write('\t'.join(map(str, outline)) + '\n')
            
            else:
                pos = int(3*j)
                outline[pos] = maxlike_dict[key][1]
                outline[pos+1] = maxlike_dict[key][0]
                outline[pos+2] = maxlike_dict[key][2]
                    
                    
                f.write('\t'.join(map(str, outline)) + '\n')
                
        f.close()
            