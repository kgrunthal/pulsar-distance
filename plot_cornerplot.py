#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 16:00:16 2023

@author: kgrunthal
"""

import numpy as np
import os, glob, json, sys, argparse
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import scipy.linalg as sl
import acor
from chainconsumer import ChainConsumer
import corner


import scipy.constants as sc

SOLAR2S = sc.G / sc.c**3 *1.98855e30
MPC2S = sc.parsec / sc.c * 1e6

###############################################################################
#   font and label settings 

plt.rcParams['font.family'] = "serif"
#plt.rcParams['font.sans-serif'] = "Times"

plt.rcParams['text.usetex']= False
plt.rcParams['xtick.labelsize'] = 10.0
plt.rcParams['ytick.labelsize'] = 10.0
plt.rcParams['axes.labelsize'] = 20.0



# names for nice labels
labels = {
          'cos_gwtheta': r'cos$\theta_\mathrm{GW}$',
          'cos_inc': r'cos$\iota$', 
          'gwphi': r'$\phi_\mathrm{GW}$',
          'log10_Mc': r'log$_{10}M_\mathrm{c}$',
          'log10_fgw': r'log$_{10}f_\mathrm{GW}$',
          'log10_h': r'log$_{10}h$',
          'phase0': r'$\Phi_0$',
          'psi': r'$\Psi$',
          'gw_log10_A': r'$\mathregular{log}_{10}A_\mathrm{GWB}$',
          'gw_gamma': r'$\gamma_\mathrm{GWB}$',
          'cgw_costheta': r'cos$\theta_\mathrm{GW}$',
          'cgw_cosinc': r'cos$i$',
          'cgw_phi': r'$\phi_\mathrm{GW}$',
          'cgw_log10_Mc': r'log$_{10}M_\mathrm{c}$',
          'cgw_log10_fgw': r'log$_{10}f_\mathrm{GW}$',
          'cgw_log10_h': r'log$_{10}h$',
          'cgw_phase0': r'$\Phi_0$',
          'cgw_psi': r'$\Psi$',
          'nmodel': r'N'
         }



def injected_values(fgw, logMc):
    Mc = (10**logMc)*SOLAR2S
    dlum = 15*MPC2S
    h = 2 * Mc**(5/3) * (np.pi*fgw)**(2/3) /dlum
    return {'cos_gwtheta': np.cos(np.pi/2),
            'gwphi': np.pi,
            'log10_h': np.log10(h),
            'log10_fgw': np.log10(fgw),
            'log10_Mc': logMc,
            'cos_inc': np.cos(np.pi/2.),
            'phase0': np.pi,
            'psi': np.pi/2.
            }


def set_up_global_options():
    parser = argparse.ArgumentParser(description='Generate sensitivity curves for PTA datasets.')
    parser.add_argument('--dir', type=str, default=None, help='Path to the chain and parameter files.')
    parser.add_argument('--chain', type=str, nargs='*', default=None,
                        help='name of the chains')
    parser.add_argument('--pars', type=int, default=None, nargs='*', help='Constrain on the parameter range')
    parser.add_argument('--names', type=str, nargs='*', default=None,
                        help='name of the chains')

    parser.add_argument('--lmc', type=float, nargs='*', default = None, help='logmc')
    parser.add_argument('--burn', type=float, default = 0.3, help='percentage of chain to cut')
    parser.add_argument('--result', type=str, default=None, help='Directory for the output')
    parser.add_argument('--parameter', type=str, default=None, help='Name of parameter file')
    parser.add_argument('--summary', action = 'store_true', help='show constrained parameter values')
    parser.add_argument('--comparison', action = 'store_true', help='show constrained parameter values')

    return parser.parse_args()




if __name__=='__main__':
    args = set_up_global_options()
    number_chains = len(args.chain)
    params = list(np.loadtxt(args.dir + args.parameter, dtype=str))
    pars = args.pars

    if args.comparison is not True:
        colors = {'8.5': [clr.to_hex('mediumpurple'), clr.to_hex('mediumorchid'), clr.to_hex('plum')],
                  '8.7': [clr.to_hex('darkslategray'), clr.to_hex('darkcyan'), clr.to_hex('darkturquoise')],
                  '9.0': [clr.to_hex('navy'), clr.to_hex('tab:blue'), clr.to_hex('skyblue')],
                  '9.5': [clr.to_hex('green'), clr.to_hex('limegreen') , clr.to_hex('lawngreen')] 
                 }

    else:
        colors = {'8.7': [clr.to_hex('firebrick'), clr.to_hex('coral')]}
    
    linestyles = ['-', '-', '-.', '--']
 
    print('\nCreating a cornerplot with {} chains \n'.format(number_chains))

    print('Setting up Chainconsumer')
    cc = ChainConsumer()
    for ii, chain in enumerate(args.chain):
        print('... on chain', ii+1)

        chain_raw = np.loadtxt(args.dir + chain)
        #np.random.shuffle(chain_raw)
        
        if args.burn != 0.:
            burn = int(args.burn*chain_raw.shape[0])
            chain_burn = chain_raw[burn:]

            corr_length, mean, sigma = acor.acor(chain_burn.T)
            chain_final = chain_burn #[::int(corr_length)]
        else:
            cl, _, _ = acor.acor(chain_raw.T)
            chain_final = chain_raw #[::int(cl/10)]
            #print('chain {} - correlation length: {}'.format(ii, cl))



        if pars is not None:
            chain_segment = np.zeros((len(chain_final[:,0]), len(pars)) )
            lbls = []
            for j, n in enumerate(pars):
                chain_segment[:, j] = chain_final[:,n]
                lbls.append(labels[params[n]])
            

            cc.add_chain(chain_segment, parameters = lbls, name=args.names[ii], color= colors[str(args.lmc[ii])][ii], linestyle=linestyles[ii])

        else:
            lbls = [labels[parname] for parname in params]
            cc.add_chain(chain_final, parameters = lbls, name=args.names[ii], color= colors[str(args.lmc[ii])][ii], linestyle=linestyles[ii])


     
    cc.configure(max_ticks = 3, tick_font_size=14, label_font_size=12, spacing=1.0,diagonal_tick_labels=True,
                 #show_contour_labels = False, contour_labels='confidence', contour_label_font_size=14,
                 #shade_gradient=[3.0], 
                 sigmas=[1,2], shade_alpha=0.6, linewidths=1.,
                 summary=args.summary, sigma2d=False, serif=True)

    out = args.dir + args.result
    
    if args.lmc is not None and np.all(np.array(args.lmc) == args.lmc[0]):
        truthvals = injected_values(22.3e-9, args.lmc[0])
        truthdict = {}
        for p in params:
            truthdict[labels[p]] = truthvals[p]
        print(truthdict)
        cc.configure_truth(color='k')
        fig = cc.plotter.plot(legend=True, filename=out, truth = truthdict)

    else:
        fig = cc.plotter.plot(legend=True, filename=out)

    fig.savefig(out, bbox_inches='tight')


    

