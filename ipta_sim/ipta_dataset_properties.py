# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 08:59:05 2024

@author: kgrun
"""


import numpy as np
import matplotlib.pyplot as plt
import json



with open('metadata/pulsar_properties.json', 'r') as tmpf:
    pulsar_dict = json.load(tmpf)
tmpf.close()

plt.figure(figsize=(8,4))


count_low = 0
low = []
count_middle = 0
middle = []
count_high = 0
high = []

for k in pulsar_dict.keys():
    tobs = (pulsar_dict[k]['end']-pulsar_dict[k]['start'])/(365.25)
    
    if (pulsar_dict[k]['dist_dm'] > 0.805) and (pulsar_dict[k]['dist_dm'] < 1.05):
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     fmt='go', ms=0.5*tobs)
        count_low +=1
        low.append(k)
        
    elif (pulsar_dict[k]['dist_dm'] > 1.36) and (pulsar_dict[k]['dist_dm'] < 1.71):
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     fmt='ro', ms=0.5*tobs)
        count_middle += 1
        middle.append(k)
        
    elif (pulsar_dict[k]['dist_dm'] > 1.9):
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     fmt='bo', ms=0.5*tobs)
        count_high += 1
        high.append(k)
        
    else:
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     marker='o',color='lightgray', ms=0.5*tobs)
    #if pulsar_dict[k]['dist_dm'] > 2:
    #    plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
    #                 fmt='bo', ms=0.5*tobs)

plt.xlabel('DM distance / kpc')
plt.ylabel('ToA error / $\mu$s')
plt.show()
print('low', count_low)
print(low)
print('middle', count_middle)
print(middle)
print('high', count_high)
print(high)




### cut out

plt.figure(figsize=(8,4))


count_low = 0

count_middle = 0
middle = []
count_high = 0
high = []
psr_all = []
psr_10 = []
psr_15 = []
psr_20 = []


for k in pulsar_dict.keys():
    tobs = (pulsar_dict[k]['end']-pulsar_dict[k]['start'])/(365.25)
    
    if (pulsar_dict[k]['dist_dm'] < 1):
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     marker = 'o', color='lightskyblue', ms=0.5*tobs)
        count_low +=1
        psr_all.append(k)
        
    elif (pulsar_dict[k]['dist_dm'] >= 1) and (pulsar_dict[k]['dist_dm'] < 1.5):
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     marker='o', color='dodgerblue',  ms=0.5*tobs)
        count_middle += 1
        psr_all.append(k)
        psr_10.append(k)
        
    elif (pulsar_dict[k]['dist_dm'] >= 1.5) and (pulsar_dict[k]['dist_dm'] < 2.0):
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     marker='o',color='mediumblue', ms=0.5*tobs)
        count_high += 1
        psr_all.append(k)
        psr_10.append(k)
        psr_15.append(k)
        
    else:
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     marker='o',color='navy', ms=0.5*tobs)
        psr_all.append(k)
        psr_10.append(k)
        psr_15.append(k)
        psr_20.append(k)
    
    
    
plt.xlabel('DM distance / kpc')
plt.ylabel('ToA error / $\mu$s')
plt.show()

print('all\n',psr_all)
print('> 1.0 kpc \n',psr_10)
print('> 1.5 kpc \n',psr_15)
print('> 2.0 kpc \n',psr_20)
