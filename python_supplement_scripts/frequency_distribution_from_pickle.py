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
import matplotlib.patches as patches

from astropy import units as u
from astropy.coordinates import SkyCoord


plt.rcParams['font.family'] = "serif"
plt.rcParams['font.sans-serif'] = "Times"

plt.rcParams['text.usetex']= False
plt.rcParams['xtick.labelsize'] = 11.0
plt.rcParams['ytick.labelsize'] = 11.0
plt.rcParams['axes.labelsize'] = 14.0



SOLAR2S = sc.G / sc.c**3 * 1.98855e30
KPC2S = sc.parsec / sc.c * 1e3
MPC2S = sc.parsec / sc.c * 1e6


'''
### FROM PICKLE ###
from enterprise_extensions import model_utils as util


def gw_frequencies_psr(psr, gwtheta, gwphi, mc, dist, fgw):
    
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

    pd = psr._pdist
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
'''  


### RAW ###

colors_85= ['mediumpurple', 'mediumorchid', 'plum']  # 0.8, 0.9, 1.0
#colors_87= ['darkred', 'red', 'orange']  # 0.8, 0.9, 1.0
colors_87= ['darkslategray', 'darkcyan', 'darkturquoise']  # 0.8, 0.9, 1.0
colors_90= ['navy', 'tab:blue', 'skyblue']  # 0.8, 0.9, 1.0
colors_95= ['green', 'limegreen' , 'lawngreen']  # 0.8, 0.9, 1.0

color_tabl = [colors_95, colors_90, colors_85]

legend_names = ['$d_\mathrm{p} = 1.0$kpc', '$d_\mathrm{p} = 1.5$kpc', '$d_\mathrm{p} = 2.0$kpc']


def gw_frequencies(ptheta, pphi, pdist, gwtheta, gwphi, mc, dist, fgw):
    
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
    ptheta = ptheta
    pphi = pphi

    # use definition from Sesana et al 2010 and Ellis et al 2012
    phat = np.array([np.sin(ptheta) * np.cos(pphi), np.sin(ptheta) * np.sin(pphi), np.cos(ptheta)])

    cosMu = -np.dot(omhat, phat)

    pd = pdist
    pd *= KPC2S  # convert from kpc to seconds

    # frequencies and phases
    omega = w0
    omega_p = w0 * (1 + fac1 * pd * (1 - cosMu)) ** (-3 / 8)
            
    return omega / np.pi, omega_p / np.pi


def generate_isotropic_distribution(Npsr):
    i = np.arange(0, Npsr, dtype=float) + 0.5
    golden_ratio = (1 + 5**0.5)/2
    costhetas = 1 - 2*i/Npsr
    thetas = np.arccos(costhetas) - np.pi/2
    phis = np.mod(2 * np.pi * i / golden_ratio, 2*np.pi)
    
    return thetas, phis


# GW source details
fgw = 22.3e-9
gw_dist = 15
gwtheta = np.pi
gwphi = np.pi


# PTA details
Npsr = 20
pthetas, pphis = generate_isotropic_distribution(Npsr)

Tobs = 10*365.25*24*3600
f_pta = 1/Tobs

bin_min = 200

fig, ax = plt.subplots(figsize=(8,5))


for ll, lmc in enumerate([9.5, 9.0, 8.5]):
    
    legend_patches = []
    for dd, pdist in enumerate([1.0, 1.5, 2.0]):
        distances = pdist*np.ones(Npsr)
        f_e, f_p = gw_frequencies(pthetas, pphis, distances, gwtheta, gwphi, 10**lmc, gw_dist, fgw)
        _,_,patch = ax.hist(f_p*1e9, bins=7,
                             color=color_tabl[ll][dd], alpha=0.8)
        legend_patches.append(patch[0])
        
        bin_min_temp = int(np.min(f_p)/f_pta)
        if bin_min_temp < bin_min:
            bin_min = bin_min_temp

    legend = ax.legend(legend_patches, legend_names,
                       title='log$_{10}M_\mathrm{c}$ = '+ '{}'.format(lmc),
                       bbox_to_anchor=(0.3+0.3*ll, -0.17),
                       frameon=False) 
    
    fig.add_artist(legend)
    
bin_max = 1 + int(f_e/f_pta)

for i in range(bin_min, bin_max):
    ax.axvline(i*f_pta*1e9, ls='--', lw = 1, color='gray', zorder=0)
    ax.text(i*f_pta*1e9, 10.5, '${}/T$'.format(i),
             fontsize=9, color='gray', backgroundcolor='w',
             ha='center')

ax.axvline(f_e*1e9, lw=3, color='red')
ax.text(f_e*1e9+0.2, 8, '$f_\mathrm{GW, Earth}$',
         color='r', rotation='vertical', fontsize=12)

# box for legend
box = patches.Rectangle((0,0), 1, 1, fill=False, edgecolor='none',
                                 visible=False)
plt.legend([box, box], ['', ''], ncol=2, borderpad=4, columnspacing=25,
           bbox_to_anchor=(0.93, -0.15))

plt.ylim(0,11.5)
plt.xlim(right=24)
plt.xlabel('$f_\mathrm{GW}$ / nHz')
plt.ylabel('number of pulsars')

#plt.savefig('frequency_distribution-isotropic.png', dpi=400, bbox_inches='tight')
plt.show()



