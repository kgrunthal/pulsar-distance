# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 18:19:39 2024

@author: kgrun
"""


import json 
import numpy as np


dictionary = json.load(open('noise.json', 'r'))

print(dictionary.keys())


k = dictionary.keys()
equad_keys = [ky for ky in k if 'EQUAD' in ky]
efac_keys = [ky for ky in k if 'EFAC' in ky]
    


psr_names = [k[:10] for k in equad_keys]

newdict = {}
for i, p in enumerate(psr_names):
    equad_new = dictionary[equad_keys[i]]*dictionary[efac_keys[i]]*1e-6
    print(p, np.log10(equad_new))
    newdict[p + '_efac'] = dictionary[efac_keys[i]]
    newdict[p + '_log10_tnequad'] = np.log10(equad_new)
    
print(newdict)


file_new = open('WN_dictionary.json', 'w')
json.dump(newdict, file_new, indent=4)
file_new.close()