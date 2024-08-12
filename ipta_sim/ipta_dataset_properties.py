# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 08:59:05 2024

@author: kgrun
"""


import numpy as np
import matplotlib.pyplot as plt
import json

from astropy import units as u
from astropy.coordinates import SkyCoord



with open('metadata/pulsar_properties.json', 'r') as tmpf:
    pulsar_dict = json.load(tmpf)
tmpf.close()

plt.figure(figsize=(8,4))


dinfo = {'low': {'npsr': 0, 'start': 60000, 'end': 0, 'tobs': 0, 'fmin':0},
         'mid': {'npsr': 0, 'start': 60000, 'end': 0, 'tobs': 0, 'fmin':0},
         'high': {'npsr': 0, 'start': 60000, 'end': 0, 'tobs': 0, 'fmin':0},
         'over_1.0': {'npsr': 0, 'start': 60000, 'end': 0, 'tobs': 0, 'fmin':0},
         'over_1.5': {'npsr': 0, 'start': 60000, 'end': 0, 'tobs': 0, 'fmin':0}
         }
count_low = 0
count_middle = 0
count_high = 0
low = []
middle = []
high = []

start_min = 60000
end_max = 0
for k in pulsar_dict.keys():
    tobs = (pulsar_dict[k]['end']-pulsar_dict[k]['start'])/(365.25)
    
    if start_min > pulsar_dict[k]['start']:
        start_min = pulsar_dict[k]['start']
    if end_max < pulsar_dict[k]['end']:
        end_max = pulsar_dict[k]['end']
    
    # low ####################################################################
    if (pulsar_dict[k]['dist_dm'] > 0.805) and (pulsar_dict[k]['dist_dm'] < 1.05):
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     fmt='go', ms=0.5*tobs)
        
        low.append(k)
        dinfo['low']['npsr'] +=1
        if dinfo['low']['start'] > pulsar_dict[k]['start']:
            dinfo['low']['start'] = pulsar_dict[k]['start']
        if dinfo['low']['end'] < pulsar_dict[k]['end']:
            dinfo['low']['end'] = pulsar_dict[k]['end']
        
        
    # middle ##################################################################
    elif (pulsar_dict[k]['dist_dm'] > 1.36) and (pulsar_dict[k]['dist_dm'] < 1.71):
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     fmt='ro', ms=0.5*tobs)
        
        middle.append(k)
        dinfo['mid']['npsr'] +=1
        if dinfo['mid']['start'] > pulsar_dict[k]['start']:
            dinfo['mid']['start'] = pulsar_dict[k]['start']
        if dinfo['mid']['end'] < pulsar_dict[k]['end']:
            dinfo['mid']['end'] = pulsar_dict[k]['end']
        
    elif (pulsar_dict[k]['dist_dm'] > 1.9):
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     fmt='bo', ms=0.5*tobs)

        high.append(k)
        dinfo['high']['npsr'] +=1
        if dinfo['high']['start'] > pulsar_dict[k]['start']:
            dinfo['high']['start'] = pulsar_dict[k]['start']
        if dinfo['high']['end'] < pulsar_dict[k]['end']:
            dinfo['high']['end'] = pulsar_dict[k]['end']
        
    else:
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     marker='o',color='lightgray', ms=0.5*tobs)


plt.xlabel('DM distance / kpc')
plt.ylabel('ToA error / $\mu$s')
plt.savefig('ipta_batches_distribution.png', bbox_inches='tight', dpi=400)
plt.show()


print(dinfo)


### cut out

plt.figure(figsize=(8,4))

psr_all = []
psr_10 = []
psr_15 = []
psr_20 = []


for k in pulsar_dict.keys():
    tobs = (pulsar_dict[k]['end']-pulsar_dict[k]['start'])/(365.25)
    
    # full data ###############################################################
    if (pulsar_dict[k]['dist_dm'] < 1):
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     marker = 'o', color='lightskyblue', ms=0.5*tobs)

        psr_all.append(k)
        
    elif (pulsar_dict[k]['dist_dm'] >= 1) and (pulsar_dict[k]['dist_dm'] < 1.5):
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     marker='o', color='dodgerblue',  ms=0.5*tobs)
        
        psr_all.append(k)
        psr_10.append(k)
        dinfo['over_1.0']['npsr'] +=1
        if dinfo['over_1.0']['start'] > pulsar_dict[k]['start']:
            dinfo['over_1.0']['start'] = pulsar_dict[k]['start']
        if dinfo['over_1.0']['end'] < pulsar_dict[k]['end']:
            dinfo['over_1.0']['end'] = pulsar_dict[k]['end']
        
    elif (pulsar_dict[k]['dist_dm'] >= 1.5) and (pulsar_dict[k]['dist_dm'] < 2.0):
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     marker='o',color='mediumblue', ms=0.5*tobs)
        
        psr_all.append(k)
        psr_10.append(k)
        psr_15.append(k)
        dinfo['over_1.0']['npsr'] +=1
        if dinfo['over_1.0']['start'] > pulsar_dict[k]['start']:
            dinfo['over_1.0']['start'] = pulsar_dict[k]['start']
        if dinfo['over_1.0']['end'] < pulsar_dict[k]['end']:
            dinfo['over_1.0']['end'] = pulsar_dict[k]['end']
        dinfo['over_1.5']['npsr'] +=1
        if dinfo['over_1.5']['start'] > pulsar_dict[k]['start']:
            dinfo['over_1.5']['start'] = pulsar_dict[k]['start']
        if dinfo['over_1.5']['end'] < pulsar_dict[k]['end']:
            dinfo['over_1.5']['end'] = pulsar_dict[k]['end']
        
    else:
        plt.errorbar(pulsar_dict[k]['dist_dm'], pulsar_dict[k]['toa_err'],
                     marker='o',color='navy', ms=0.5*tobs)
        psr_all.append(k)
        psr_10.append(k)
        psr_15.append(k)
        psr_20.append(k)
    
    
    
plt.xlabel('DM distance / kpc')
plt.ylabel('ToA error / $\mu$s')
plt.savefig('ipta_dropout_distribution.png', bbox_inches='tight', dpi=400)
plt.show()



### gather stats ##############################################################

for key in dinfo.keys():
    tobs = dinfo[key]['end']-dinfo[key]['start']
    fmin = 1/(tobs*24*3600)
    dinfo[key]['tobs'] = tobs/365.25
    dinfo[key]['fmin'] = fmin
    
'''
print('all\n',psr_all)
print('> 1.0 kpc \n',psr_10)
print('> 1.5 kpc \n',psr_15)
print('> 2.0 kpc \n',psr_20)
'''

f_min = 1/((end_max-start_min)*24*3600)
print('FULL DATASET:')
print('    fundamental frequency = ', f_min)
print('    bin of CGW', 22e-9/f_min)
print()

print('low dataset')
print('    number of pulsars: {}'.format(dinfo['low']['npsr']))
print('    Tobs: {} yr'.format(dinfo['low']['tobs']))
print('    fmin: {} Hz'.format(dinfo['low']['fmin']))
print()

print('mid dataset')
print('    number of pulsars: {}'.format(dinfo['mid']['npsr']))
print('    Tobs: {} yr'.format(dinfo['mid']['tobs']))
print('    fmin: {} Hz'.format(dinfo['mid']['fmin']))
print()

print('high dataset')
print('    number of pulsars: {}'.format(dinfo['high']['npsr']))
print('    Tobs: {} yr'.format(dinfo['high']['tobs']))
print('    fmin: {} Hz'.format(dinfo['high']['fmin']))
print()

print('> 1.0 kpc dataset')
print('    number of pulsars: {}'.format(dinfo['over_1.0']['npsr']))
print('    Tobs: {} yr'.format(dinfo['over_1.0']['tobs']))
print('    fmin: {} Hz'.format(dinfo['over_1.0']['fmin']))

print('> 1.5 dataset')
print('    number of pulsars: {}'.format(dinfo['over_1.5']['npsr']))
print('    Tobs: {} yr'.format(dinfo['over_1.5']['tobs']))
print('    fmin: {} Hz'.format(dinfo['over_1.5']['fmin']))




#### pulsar positions #########################################################

plt.figure(figsize=(7,9))
plt.subplot(projection="mollweide")

for psr in low:
    name = 'PSR '+psr
    pos = SkyCoord.from_name(name, parse=True)
    ra, dec = pos.ra.wrap_at(180*u.degree).radian, pos.dec.radian
    
    plt.errorbar(ra, dec, fmt='go')
    
    print(name, end=' \r')


for psr in middle:
    name = 'PSR '+psr
    pos = SkyCoord.from_name(name, parse=True)
    ra, dec = pos.ra.wrap_at(180*u.degree).radian, pos.dec.radian
    
    plt.errorbar(ra, dec, fmt='ro')
    
    print(name, end=' \r')
    
    
for psr in high:
    name = 'PSR '+psr
    pos = SkyCoord.from_name(name, parse=True)
    ra, dec = pos.ra.wrap_at(180*u.degree).radian, pos.dec.radian
    
    plt.errorbar(ra, dec, fmt='bo')
    
    print(name, end=' \r')
    
plt.grid(which='both')
plt.savefig('ipta_batches_skydistribution.png', bbox_inches='tight', dpi=400)
plt.show()






plt.figure(figsize=(7,9))
plt.subplot(projection="mollweide")

for psr in psr_all:
    name = 'PSR '+psr
    pos = SkyCoord.from_name(name, parse=True)
    ra, dec = pos.ra.wrap_at(180*u.degree).radian, pos.dec.radian
    
    plt.errorbar(ra, dec, marker='*', color='lightskyblue')
    
    print(name, end=' \r')
    
for psr in psr_10:
    name = 'PSR '+psr
    pos = SkyCoord.from_name(name, parse=True)
    ra, dec = pos.ra.wrap_at(180*u.degree).radian, pos.dec.radian
    
    plt.errorbar(ra, dec, marker='*', color='dodgerblue')
    
    print(name, end=' \r')


for psr in psr_15:
    name = 'PSR '+psr
    pos = SkyCoord.from_name(name, parse=True)
    ra, dec = pos.ra.wrap_at(180*u.degree).radian, pos.dec.radian
    
    plt.errorbar(ra, dec, marker='*', color='mediumblue')
    
    print(name, end=' \r')
    
    
for psr in high:
    name = 'PSR '+psr
    pos = SkyCoord.from_name(name, parse=True)
    ra, dec = pos.ra.wrap_at(180*u.degree).radian, pos.dec.radian
    
    plt.errorbar(ra, dec, marker='*', color='navy')
    
    print(name, end=' \r')
    
plt.grid(which='both')
plt.savefig('ipta_dropout_skydistribution.png', bbox_inches='tight', dpi=400)
plt.show()





