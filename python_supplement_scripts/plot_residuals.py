#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 14:25:56 2024

@author: kgrunthal
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np


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


#basepath = '/u/kgrunthal/HD/MCMCout_RNCGW9.0_zeta1.0/'
basepath = '/u/kgrunthal/HD/MCMCout_GWB_test/'
psrs = pickle.load(open(basepath + 'psrs.pkl', 'rb'))


for p in psrs:
    print(p.name)
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.errorbar(p.toas, p.residuals, yerr=p.toaerrs,capsize=2, ls='', fmt='kx')
    ax2.plot(p.toas, p.residuals, ls='')
    #print(ax1.get_yticks())
    ax2.set_yticks(ax1.get_yticks(), labels= [np.round(x, 2) for x in ax1.get_yticks()*F0[p.name]*1e3] )
    
    plt.title(p.name)
    ax1.set_xlabel('MJD')
    ax1.set_ylabel('residual / s')
    ax2.set_ylabel('$\Delta\Phi$ / 1e-3')
    plt.savefig(basepath + p.name + '.png', bbox_inches='tight')
    
