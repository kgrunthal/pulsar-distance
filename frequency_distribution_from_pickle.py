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

    pd = pdist

    # convert units
    pd *= KPC2S  # convert from kpc to seconds

    # frequencies and phases
    omega = w0
    omega_p = w0 * (1 + fac1 * pd * (1 - cosMu)) ** (-3 / 8)
            
    return omega, omega_p




def set_up_global_options():
    parser = argparse.ArgumentParser(description='Generate sensitivity curves for PTA datasets.')
    parser.add_argument('--folder', type=str, default=None, help='Path to the MCMCfolder.')
    parser.add_argument('--lmc', type=float, nargs='*', default=9.5)
    parser.add_argument('--fgw', type=float, nargs='*', default=22.3)

    return parser.parse_args()




args = set_up_global_options()

with open(args.folder+'psrs.pkl', 'rb') as psrfile:
    PSRs =  pickle.load(psrfile)
psrfile.close()



gwtheta = np.pi
gwphi = np.pi
gw_dist = 15   #Mpc


dtype = [('name', str), ('freq', float)]
tuples = []
for ii, psr in enumerate(PSRs):
    omega, omega_p = gw_frequencies(psr, gwtheta, gwphi, 10**args.lmc, gw_dist, args.fgw*10**(-9))
    tuples.append((psr.name, omega_p))
    
values = np.array(tuples, dtype)

sorted_values = np.sort(values, order = 'freq')
print(sorted_values)
    

    
    
