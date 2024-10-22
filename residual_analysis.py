#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 09:00:32 2024

@author: kgrunthal
"""

import numpy as np
import pickle 
import argparse
import scipy.fft as scfft
import matplotlib.pyplot as plt


def set_up_global_options():
    parser = argparse.ArgumentParser(description='CGW analysis for MeerKAT PTA')
    
    parser.add_argument('--psrpickle', type=str, default=None,
                        help='pickle file with PSRs object')
    
    parser.add_argument('--outdir', type=str, default=None,
                        help='Directory for the MCMC output folder')
    
    return parser.parse_args()

F0 = {'J0024-2029': 238.6718614640,
      'J0114+7148': 276.9603731761,
      'J0234+0837': 372.8217959796,
      'J0354-4032': 208.5405026916,
      'J0444+4032': 345.8937441534,
      'J0604-0837': 483.6205702139,
      'J0724-7148': 236.6732752649,
      'J0814+2029': 374.2537559448,
      'J0934-2644': 581.3430326118,
      'J1024+5812': 587.1417547818,
      'J1144+0251': 675.6059855957,
      'J1304-4835': 999.0770504809,
      'J1354+3322': 487.4976005103,
      'J1514-1428': 775.8871495052,
      'J1724+1428': 535.4625206295,
      'J1844-3322': 235.9218726028,
      'J1934+4835': 354.7123046593,
      'J2054-0251': 311.1677553594,
      'J2214-5812': 992.8428043529,
      'J2304+2644': 282.4822855868
      }



def main():
    args = set_up_global_options()
    print('LOAD DATA')
    print('... loading pulsars from pickle \n')
    with open(args.psrpickle, 'rb') as psrpickle:
        PSRs = pickle.load(psrpickle)
    psrpickle.close()
    
    for psr in PSRs:
        print('on PSR ', psr.name)
        print('... fourier analysis of residuals')
        res = psr.residuals
        N = len(res)
        dt  = (14.*24.*3600.)
        xf = scfft.fftfreq(N, dt)
        yf = scfft.fft(res)
        print('... plotting residuals')
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize = (14,5) )
        
        ax0_2 = axs[0].twinx()
        axs[0].errorbar(psr.toas, res, yerr=psr.toaerrs, capsize=2, ls='', fmt='kx')
        ax0_2.plot(psr.toas, res, ls='')
        
        axs[1].errorbar(np.log10(xf[1:N//2]), 2.0/N * np.abs(yf[1:N//2]), ls='', fmt='kx')
        
        plt.suptitle(psr.name)
        
        axs[0].set_xlabel('MJD')
        axs[0].set_ylabel('residual / s')
        ax0_2.set_yticks(axs[0].get_yticks(), labels= [np.round(x, 2) for x in axs[0].get_yticks()*F0[psr.name]*1e3] )        
        ax0_2.set_ylabel('$\Delta\Phi$ / 1e-3')

        axs[1].set_xlabel('frequency')
        plt.savefig(args.outdir + psr.name + '.png', bbox_inches='tight')
    
    return 0


if __name__ == '__main__':
    main()
        
