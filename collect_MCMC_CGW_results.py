# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 13:41:10 2024

@author: kgrun
"""

import json
import glob, os

import numpy as np

def write_outline(openfile, maxdict, nparameter, idx):
    outline = np.zeros(3*nparameter)
    if i == 1:
        headline =''
        for jj, key in enumerate(maxlike_dict.keys()):
            if jj == 0:
                headline += '{}\t{}_ll\t{}_ul'.format(key, key, key)
            else:
                headline += '\t{}\t{}_ll\t{}_ul'.format(key, key, key)

            pos = int(3*jj)
            outline[pos] = maxlike_dict[key][1]
            outline[pos+1] = maxlike_dict[key][0]
            outline[pos+2] = maxlike_dict[key][2]

        openfile.write(headline + '\n')
        openfile.write('\t'.join(map(str, outline)) + '\n')

    else:
        for jj, key in enumerate(maxlike_dict.keys()):
            pos = int(3*jj)
            outline[pos] = maxlike_dict[key][1]
            outline[pos+1] = maxlike_dict[key][0]
            outline[pos+2] = maxlike_dict[key][2]

        openfile.write('\t'.join(map(str, outline)) + '\n')

    return None
 
basepath='/u/kgrunthal/HD/'
outdir = '/out/CGWparameter_recovery/ipta/RA18hDEC-45deg/high/'

start = 1
end = 50

#for lmc in [8.6]:
for lmc in [8.5, 9.0, 9.5]:
    for pd in [1.0]:
        print(lmc, ' ', pd, ':')
        
        # create outputfiles
        f_noPT = open(basepath + outdir + '/parameters_noPT_lmc{}_pd{}.txt'.format(lmc, pd), 'a')
        #f_PT_pdist = open(basepath+'/out/CGWparameter_recovery/isotropic/parameters_PT_pdist_lmc{}_pd{}.txt'.format(lmc, pd), 'a')
        #f_PT_pphase = open(basepath+'/out/CGWparameter_recovery/isotropic/parameters_PT_pphase_lmc{}_pd{}.txt'.format(lmc, pd), 'a')
        #f_PT_pdist_pphase = open(basepath+'/out/CGWparameter_recovery/isotropic/parameters_PT_lmc{}_pd{}.txt'.format(lmc, pd), 'a')

        #iterate over realisations
        for i in range(start,end+1):
            #path = basepath+'/MCMCout_IPTA_CGW{}_pd{}_{}'.format(lmc, pd, i)
            path = basepath+'/MCMCout_IPTA_CGW_RA18h-45deg_{}_high_{}'.format(lmc, i)
            #path = basepath+'/MCMCout_IPTAhigh_CGW{}_mid_{}'.format(lmc, i)
            print('... on {}'.format(path))
            
            ## no pulsar term #################################################
            if os.path.isfile(path + '/CGWsearch_noPT/maxlike.json'):
                with open(path + '/CGWsearch_noPT/maxlike.json', 'r') as noPT_file:
                    maxlike_dict = json.load(noPT_file)
                noPT_file.close()
                write_outline(f_noPT, maxlike_dict, 8, i)

            
            '''
            ## pulsar term - pdist #############################################

            if os.path.isfile(path + '/CGWsearch_PT_pdist/maxlike.json'):
                with open(path + '/CGWsearch_PT_pdist/maxlike.json', 'r') as PT_file:
                    maxlike_dict = json.load(PT_file)
                PT_file.close()
                write_outline(f_PT_pdist, maxlike_dict, 28, i)

             
            ## pulsar term - pphase ############################################

            if os.path.isfile(path + '/CGWsearch_PT_pphase/maxlike.json'):
                with open(path + '/CGWsearch_PT_pphase/maxlike.json', 'r') as PT_file:
                    maxlike_dict = json.load(PT_file)
                PT_file.close()
                write_outline(f_PT_pphase, maxlike_dict, 28, i)

            

            ## pulsar term - pdist-pphase  ######################################

            if os.path.isfile(path + '/CGWsearch_PT/maxlike.json'):
                with open(path + '/CGWsearch_PT/maxlike.json', 'r') as PT_file:
                    maxlike_dict = json.load(PT_file)
                PT_file.close()
                write_outline(f_PT_pdist_pphase, maxlike_dict, 48, i)

            '''

        f_noPT.close()
        #f_PT_pdist.close()
        #f_PT_pphase.close()
        #f_PT_pdist_pphase.close()
