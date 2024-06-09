# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 16:48:15 2024

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

import ephem

import libstempo.toasim as LT
import libstempo.plot as LP

from enterprise.signals import signal_base, gp_signals, white_signals, utils, gp_priors
from enterprise.pulsar import Pulsar
import enterprise.signals.parameter as parameter

from enterprise_extensions import model_utils, blocks, deterministic
from enterprise_extensions.frequentist import optimal_statistic as opt_stat
from enterprise_extensions import sampler

import corner
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc




def set_up_global_options():
    parser = argparse.ArgumentParser(description='Generate sensitivity curves for PTA datasets.')
    parser.add_argument('--basedir', type=str, default=None, help='Path to the parfiles.')
    parser.add_argument('--outdir', type=str, default=None, help='Path to the parfiles.')
    
    parser.add_argument('--ptamodel', type=str, default='CGW', help='Which signals to include in the PTA model')
    parser.add_argument('--lmc', type=float, nargs='*', default=9.5)
    parser.add_argument('--fgw', type=float, nargs='*', default=22.3)
    parser.add_argument('--ncgw', type=int, default=1)

    parser.add_argument('--psrTerm', action="store_true", default=False, help='Simulate CGW with or without pulsar term')
    parser.add_argument('--zeta', type=float, default=1.0)
    
    parser.add_argument('--analysis', action="store_true", default=False, help='Save cornerplot from the MCMC run')
    
    parser.add_argument('--N', type=int, default=None, help='number of MCMC draws')
    return parser.parse_args()



def PTA_model(psrs, ptamodel, psrTerm=False):
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
            
        elif m=='CGW':
            cos_gwtheta = parameter.Uniform(-1, 1)                  # position of source
            gwphi = parameter.Uniform(0, 2*np.pi)                   # position of source
            log10_mc = parameter.Uniform(6.5, 10.)                  # chirp mass
            log10_h = parameter.Uniform(-16, -11)                   # strain amplitude
            log10_fgw = parameter.Uniform(-8., -7.)                 # gw frequency
            phase0 = parameter.Uniform(0, 2*np.pi)                  # gw phase
            psi = parameter.Uniform(0, np.pi)                       # gw polarization 
            cos_inc = parameter.Uniform(-1, 1)                      # inclination of binary with respect to Earth 
            
            cw_wf = deterministic.cw_delay(cos_gwtheta=cos_gwtheta, gwphi=gwphi, log10_mc=log10_mc, 
                                           log10_h=log10_h, log10_fgw=log10_fgw, phase0=phase0, 
                                           psi=psi, cos_inc=cos_inc)
            CGW = deterministic.CWSignal(cw_wf, ecc=False, psrTerm=False)
            s.append(CGW)
            outstring += ' CGW '
 
        else:
            print('Unsupported PTA model component')
            break
    
    signal = s[0]
    for sl in s[1:]:
        signal+=sl
    
    print(outstring, flush=True)






def produce_output(outdir = ''):
    chain = np.loadtxt(outdir + '/chain_1.txt')
    
    params = []
    with open(outdir + '/pars.txt', 'r') as file:
        for line in file.readlines():
            params.append(line.strip())
    


    burn = int(0.3*chain.shape[0])   # burn = beginning of chain, not helpful!
    outdict={}
    
    for i, p in enumerate(params):    
        # sigma
        parameter_estimate_list = corner.quantile(chain[burn:,i], [0.16, 0.5, 0.84])
        # replace mean with maximum likelihood
        n, bins, patches = plt.hist(chain[burn:,i], bins=100)
        plt.clf()
        max_lh_pos = np.argmax(n)
        max_lh = bins[max_lh_pos]
        parameter_estimate_list[1] = max_lh
        
        outdict[p]=parameter_estimate_list
        
    
    corner.corner(chain[burn:,:-4], labels=params[:],
                  bins =30,
                  plot_datapoints=False, plot_density=True, 
                  plot_contours=False,fill_contours=False,
                  show_titles = True, use_math_text=True, quantiles=[0.16, 0.5, 0.84], verbose=True)
    plt.savefig('{}/cornerplot.png'.format(outdir), bbox_inches='tight')
    plt.clf()
        
    with open(outdir+"/maxlike.json", "w") as f:
        json.dump(outdict, f, indent=4)   
    f.close()
    
    return chain[burn:]




def main():
    args = set_up_global_options
    
    outD = args.outdir
    
    with open(args.basedir + '/psrs.pkl', 'rb') as psrpickle:
        ePSRs = pickle.load(psrpickle)
        print('loaded pulsars from pickle')
    psrpickle.close()
   
    
    PTA = PTA_model(ePSRs, args.ptamodel, psrTerm=args.psrTerm)
    
    
    x0 = np.hstack([p.sample() for p in PTA.params])
    iteration = 0
    print('... looking for a suitable initial sample, iteration {}'.format(iteration), end = '\r', flush=True)
    while PTA.get_lnlikelihood(x0) < 0. or PTA.get_lnprior(x0) == float(-np.inf):
        iteration += 1
        x0 = np.hstack([p.sample() for p in PTA.params])
        print('\r ... looking for a suitable initial sample, iteration {}'.format(iteration), end = '\r',flush=True)
    print(' \n ... found a sample' , flush = True)

    
    smpl = sampler.setup_sampler(PTA, outdir=outD)
    smpl.sample(x0, args.Nsample, SCAMweight=60, AMweight=0, DEweight=30, burn=int(0.3*args.Nsample))
    
    if args.analysis == True:
        _ = produce_output(outdir = outD)
        return None
    
    else:
        return None
    

if __name__ == '__main__':
    main()
        