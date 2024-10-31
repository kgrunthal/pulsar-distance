#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 11:35:11 2024

@author: kgrunthal
"""

import numpy as np
import os, glob, json, sys, argparse
from pathlib import Path
import matplotlib
#matplotlib.use('GTK3Agg')
import matplotlib.pyplot as plt
import scipy.linalg as sl
import acor
#from chainconsumer import ChainConsumer
import corner








def set_up_global_options():
    parser = argparse.ArgumentParser(description='Generate sensitivity curves for PTA datasets.')
    parser.add_argument('--basestring', type=str, default=None, help='Path to the chain and parameter files.')
    parser.add_argument('--path', type=str, default=None, help='Path to the chain and parameter files.')
    parser.add_argument('--outfile', type=str, default=None,  help='name of the chain')
    return parser.parse_args()




if __name__=='__main__':
    args = set_up_global_options()
    analysis_folder = '/CGWsearch_noPT/' 
    MCMCfolder_list = sorted(glob.glob(args.path + args.basestring + '*'))
    print(MCMCfolder_list)
    MCMCfolder = []
    for folder in MCMCfolder_list:
        chainfile = Path(folder + analysis_folder + '/chain_1.txt')
        if chainfile.exists():
            MCMCfolder.append(folder)
        else:
            print('Skip', folder)

    # initialise full chain with first folder
    chain_raw = np.loadtxt(MCMCfolder[0] + analysis_folder + '/chain_1.txt')
    
    burn = int(0.3*chain_raw.shape[0])
    chain = chain_raw[burn:]

    corr_length, mean, sigma = acor.acor(chain.T)
    chain_thinned = chain[::int(corr_length)]
    print(MCMCfolder[0],' ... correlation length: ', corr_length, ' ... thinned chain length: ', len(chain_thinned))
    
    chain_full = chain_thinned.copy()
    
    for folder in MCMCfolder[1:]:
        chain_raw = np.loadtxt(folder + analysis_folder + '/chain_1.txt')
        
        burn = int(0.3*chain_raw.shape[0])
        chain = chain_raw[burn:]

        corr_length, mean, sigma = acor.acor(chain.T)
        if int(corr_length) == 0:
            slc = 1
        else:
            slc = int(corr_length)
        chain_thinned = chain[::int(slc)]
        print(folder,' ... correlation length: ', corr_length, ' ... thinned chain length: ', len(chain_thinned))
        chain_full = np.append(chain_full, chain_thinned, axis = 0)
    
    np.savetxt(args.outfile, chain_full)
    

        
    
    
    
