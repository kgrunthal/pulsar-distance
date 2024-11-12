# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 13:11:56 2024

@author: kgrun
"""


import numpy as np
import json
import matplotlib.pyplot as plt


with open('distances.json', 'r') as dist_file:
    dist_dict = json.load(dist_file)
dist_file.close()

names = np.loadtxt('ska_far.txt', dtype=str)

with open('../ipta_sim/metadata/pulsar_properties.json', 'r') as tmpf:
    ipta_dict = json.load(tmpf)
tmpf.close()

with open('../ipta_sim/WN_dictionary.json', 'r') as tmpfile:
    WN_ipta = json.load(tmpfile)
tmpfile.close()

EFACs = []
EQUADs = []
for key in ipta_dict.keys():
    if ipta_dict[key]['pta'] == 'mpta':
        EFACs.append(WN_ipta[key+'_efac'])
        EQUADs.append(WN_ipta[key+'_log10_tnequad'])
    
np.random.seed(42)
cads = np.random.uniform(12, 17, 114)
efacs = np.random.uniform(np.min(np.array(EFACs)), np.max(np.array(EFACs)), 114)
equads = np.random.uniform(np.min(np.array(EQUADs)), np.max(np.array(EQUADs)), 114)

#toa_errs = np.random.uniform(4, 6, 114)
toa_errs = np.random.uniform(0.1, 4, 114)

temptoaerr = np.array([ipta_dict[key]['toa_err'] for key in ipta_dict.keys()]) + np.random.uniform(-0.1,0.1)
toa_errs = np.array([np.max([0.1, x]) for x in temptoaerr])
print(toa_errs)

SKA_dict = {}
WN_dict = {}

for n, name in enumerate(names):

    SKA_dict[name] = {"pta": "ska",
                      "toa_err": np.round(toa_errs[n], 2),
                      "start": 62137.0,
                      "end": 64328.0,
                      "cad": np.round(cads[n], 2),
                      "dist_dm": np.round(dist_dict[name], 3),
                      "dist": np.round(dist_dict[name], 3)}
    WN_dict[name+'_efac'] = efacs[n]
    WN_dict[name+'_log10_tnequad'] = equads[n]
    

with open('SKA_properties.json', 'w') as skaprop:
    json.dump(SKA_dict, skaprop, indent=4)
skaprop.close()

'''
with open('SKA_WN_dictionary.json', 'w') as skawn:
    json.dump(WN_dict, skawn, indent=4)
skawn.close()
'''




    
