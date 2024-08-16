# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 09:37:49 2024

@author: kgrun
"""



import argparse
import numpy as np
import pickle
import json
import glob
import os
import scipy.constants as sc
#import scipy.signal as signal
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord

from enterprise_extensions import model_utils as util

SOLAR2S = sc.G / sc.c**3 * 1.98855e30
KPC2S = sc.parsec / sc.c * 1e3
MPC2S = sc.parsec / sc.c * 1e6


def gw_frequencies(psr, gwtheta, gwphi, mc, dist, fgw):
    
    # convert units
    mc *= SOLAR2S  # convert from solar masses to seconds
    dist *= MPC2S  # convert from Mpc to seconds

    # define initial orbital frequency
    w0 = np.pi * fgw

    # define variable for later use
    cosgwtheta, cosgwphi = np.cos(gwtheta), np.cos(gwphi)
    singwtheta, singwphi = np.sin(gwtheta), np.sin(gwphi)

    # unit vectors to GW source
    omhat = np.array([-singwtheta * cosgwphi, -singwtheta * singwphi, -cosgwtheta])

    # various factors invloving GW parameters
    fac1 = 256 / 5 * mc ** (5 / 3) * w0 ** (8 / 3)


    # pulsar location
    ptheta = psr.theta
    pphi = psr.phi
    

    # use definition from Sesana et al 2010 and Ellis et al 2012
    phat = np.array([np.sin(ptheta) * np.cos(pphi), np.sin(ptheta) * np.sin(pphi), np.cos(ptheta)])

    cosMu = -np.dot(omhat, phat)

    pd = psr._pdist[0]
    # convert units
    pd *= KPC2S  # convert from kpc to seconds

    # frequencies and phases
    omega = w0
    omega_p = w0 * (1 + fac1 * pd * (1 - cosMu)) ** (-3 / 8)
            
    return omega / np.pi, omega_p / np.pi




def set_up_global_options():
    parser = argparse.ArgumentParser(description='Generate sensitivity curves for PTA datasets.')
    parser.add_argument('--folder', type=str, default=None, help='Path to the MCMCfolder.')
    parser.add_argument('--lmc', type=float, default=9.5)
    parser.add_argument('--fgw', type=float, default=22.3)

    return parser.parse_args()




args = set_up_global_options()

print(args.folder.split('/')[-1])
with open(args.folder+'/psrs.pkl', 'rb') as psrfile:
    PSRs =  pickle.load(psrfile)
psrfile.close()


with open(args.folder + '/distances.json') as distfile:
    pdistances = json.load(distfile)
distfile.close()



gwtheta = np.pi
gwphi = np.pi
gw_dist = 15   #Mpc

dtype = [('name', 'U11'), ('freq', float)]
tuples = []

print('calculating frequencies')
for ii, psr in enumerate(PSRs):
    psr._pdist = (pdistances[psr.name], 0.2)
    f_e, f_p = gw_frequencies(psr, gwtheta, gwphi, 10**args.lmc, gw_dist, args.fgw*10**(-9))
    tuples.append((psr.name, f_p))

values = np.array(tuples, dtype)
sorted_values = np.flipud(np.sort(values, order = 'freq'))

print('getting PTA information')

Tobs = util.get_tspan(PSRs)
f_pta = 1/Tobs

bin_min = np.round(np.min(sorted_values['freq'])/f_pta)
bin_max = np.round(f_e/f_pta)
print(bin_min, bin_max)

for jj in range(len(PSRs)):
    plt.errorbar(jj, sorted_values['freq'][jj]*1e9, fmt='ko', ls='')

plt.axhline(f_e*1e9, ls='-', color='darkgray')
plt.text(len(PSRs), f_e*1e9+0.2, '$f_e$', fontsize=7, color='darkgray')
for ff in range(int(bin_min), int(bin_max)):
    plt.axhline(ff*f_pta*1e9, ls=':', color='gray')
    plt.text(len(PSRs), ff*f_pta*1e9+0.2, '${}/T$'.format(ff), fontsize=7, color='lightgray')

plt.xticks(range(len(PSRs)), sorted_values['name'])
plt.tick_params(axis='x', labelrotation=90)
plt.ylabel('f / nHz')

plt.xlim(-1, len(PSRs)+2)
plt.savefig('./frequency_distribution/' + args.folder.split('/')[-1] + '.png', dpi=400, bbox_inches='tight')
    

    
    
