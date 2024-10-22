#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 12:09:11 2023

@author: kgrunthal
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

import ephem

import libstempo.toasim as LT
import libstempo.plot as LP

from enterprise.signals import signal_base, gp_signals, white_signals, utils, gp_priors
from enterprise.pulsar import Pulsar
import enterprise.signals.parameter as parameter

from enterprise_extensions import model_utils, blocks, deterministic
from enterprise_extensions.frequentist import optimal_statistic as opt_stat
from enterprise_extensions import sampler as sp

import corner
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc



SOLAR2S = sc.G / sc.c**3 * 1.98855e30
KPC2S = sc.parsec / sc.c * 1e3
MPC2S = sc.parsec / sc.c * 1e6




def S(A, gamma, f):
    fc = 1/(365.25*86000)
    A_c = (A**2)/(12*np.pi**2.)
    return A_c*((f/fc)**(-1*gamma))*(fc**(-3.))

def singlebin(amp, gamma, f_window, f):
    spectrum = 1e-50*np.ones(len(f))
    window = [ii for ii in range(len(f)) if f[ii] > f_window[0] and f[ii] < f_window[1]]
    spectrum[window] = S(amp, gamma, f[window])
    return spectrum


def PTA_model(psrs, obtimes, ptamodel, mode=None):
    #find the maximum time span to set GW frequency sampling
    Tspan = model_utils.get_tspan(psrs)
    
    s = []
    model = ptamodel.split(',')
    outstring = 'PTA model with '
    
    for m in model:
        if m =='TM':
            s.append(gp_signals.TimingModel())   # First we add the timing model
            outstring += 'TM '

        elif m=='WN':
            s.append(white_signals.MeasurementNoise(efac=1.))   #Timing noise
            outstring += 'WN '
    
        elif m=='RN':
            #log10_A = parameter.Uniform(-15., -12.)
            #gamma = parameter.Uniform(1.5, 2.5)
            #pl = utils.powerlaw(log10_A=log10_A, gamma=gamma)
            #s.append(gp_signals.FourierBasisGP(spectrum=pl, Tspan=Tspan, components=20) )
            s.append(blocks.red_noise_block(components = 20))
            outstring += 'RN '
        
        elif m=='CRN_fs':
            s.append(blocks.common_red_noise_block(psd='spectrum', prior='log-uniform', name='gw', Tspan = Tspan, components=20))
            outstring += 'CRN (spectrum) '

        elif m=='CRN_pl':
            s.append(blocks.common_red_noise_block(psd='powerlaw', prior='log-uniform', Tspan = Tspan, name='gw', components = 20))
            outstring += 'CRN (powerlaw) '
 
        else:
            print('Unsupported PTA model component')
            break
    
    signal = s[0]
    for sl in s[1:]:
        signal+=sl
    
    print(outstring, flush=True)
   

    # We set up the PTA object using the signal we defined above and the pulsars
    pta = signal_base.PTA([signal(p) for p in psrs])
    
    return pta




def set_up_global_options():
    parser = argparse.ArgumentParser(description='Generate sensitivity curves for PTA datasets.')
    parser.add_argument('--psrpickle', type=str, default=None, help='psrpickle')
    return parser.parse_args()






    
    
def main():
    args = set_up_global_options()
    ePSRs = pickle.load(open(args.psrpickle, 'rb'))
    
    obstimes = np.arange(50000,53652,14)
    frequencies = np.linspace(1/(3625*86400), 20/(3625*86400), 20)

    pta_pl = PTA_model(ePSRs, obstimes, "TM,WN,CRN_pl")
    pta_fs = PTA_model(ePSRs, obstimes, "TM,WN,CRN_fs")
    
    parameter_powerlaw = {"gw_gamma": 13./3.,
                          "gw_log10_A": np.log10(2e-15)}

    parameter_fs = {"gw_log10_rho_0": -5.782842379847864,
                    "gw_log10_rho_1": -6.085640453052354,
                    "gw_log10_rho_2": -6.267303605106803,
                    "gw_log10_rho_3": -6.395972028095045,
                    "gw_log10_rho_4": -6.4918242421336565,
                    "gw_log10_rho_5": -6.569342307987413,
                    "gw_log10_rho_6": -6.642622301866741,
                    "gw_log10_rho_7": -6.702933978477164,
                    "gw_log10_rho_8": -6.748473481810571,
                    "gw_log10_rho_9": -6.796641381996851,
                    "gw_log10_rho_10": -6.838228735891385,
                    "gw_log10_rho_11": -6.868593584543247,
                    "gw_log10_rho_12": -6.911613536514633,
                    "gw_log10_rho_13": -6.933384795378554,
                    "gw_log10_rho_14": -6.975868631760214,
                    "gw_log10_rho_15": -7.006216284959804,
                    "gw_log10_rho_16": -7.0268233372362365,
                    "gw_log10_rho_17": -7.054332986750679,
                    "gw_log10_rho_18": -7.076251467122198,
                    "gw_log10_rho_19": -7.091965005316915
                    }
    

    free_spectrum = np.log10(S(2e-15, 13/3, frequencies))
    #free_spectrum = singlebin(2e-15, 13/3, [3./(3625*86400), 4./(3625*86400)], frequencies)

    for ii,key in enumerate(parameter_fs.keys()):
       parameter_fs[key] = free_spectrum[ii]

    pardict_fs = pta_fs.map_params(list(parameter_fs.values()))
    pardict_fs.update(parameter_fs)


    print('OS stats powerlaw')
    
    ostat_pl = opt_stat.OptimalStatistic(ePSRs, pta=pta_pl, orf='hd')
    _, _, _, OS_pl, OSsig_pl = ostat_pl.compute_os(params=parameter_powerlaw, psd='powerlaw')

    print('OS = {} +/- {} \t S/N = {}'.format(OS_pl, OSsig_pl, OS_pl/OSsig_pl))
    print()

    
    print('OS stats free spectrum')
    ostat_fs = opt_stat.OptimalStatistic(ePSRs, pta=pta_fs, orf='hd')
    OS_fs, OSsig_fs = np.zeros(20), np.zeros(20)
    for ii, f in enumerate(frequencies):
        _, _, _, OS_fs[ii], OSsig_fs[ii] = ostat_fs.compute_os(params=pardict_fs, psd='spectrum', fgw=f)
        print('OS = {} +/- {} \t S/N = {}'.format(OS_fs[ii], OSsig_fs[ii], OS_fs[ii]/OSsig_fs[ii]))
   
    return None

    

if __name__ == '__main__':
    main()

    

