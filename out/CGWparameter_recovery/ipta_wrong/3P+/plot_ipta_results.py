# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 09:59:38 2024

@author: kgrun
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import scipy.constants as sc

SOLAR2S = sc.G / sc.c**3 *1.98855e30
MPC2S = sc.parsec / sc.c * 1e6

def injected_values(fgw, logMc):
    Mc = (10**logMc)*SOLAR2S
    dlum = 15*MPC2S
    h = 2 * Mc**(5/3) * (np.pi*fgw)**(2/3) /dlum
    return {'cos_gwtheta': np.cos(np.pi/2),
            'gwphi': np.pi,
            'log10_h': np.log10(h),
            'log10_fgw': np.log10(fgw),
            'log10_Mc': logMc,
            'cos_inc': np.cos(np.pi),
            'phase0': np.pi,
            'psi': np.pi
            }

plt.rcParams['font.family'] = "serif"
plt.rcParams['font.sans-serif'] = "Times"

plt.rcParams['text.usetex']= False
plt.rcParams['xtick.labelsize'] = 11.0
plt.rcParams['ytick.labelsize'] = 11.0
plt.rcParams['axes.labelsize'] = 14.0




colors_85= ['plum', 'mediumorchid', 'mediumpurple', 'pink']  # 0.8, 0.9, 1.0
colors_90= ['skyblue', 'tab:blue', 'navy']  # 0.8, 0.9, 1.0
colors_95= ['lawngreen', 'limegreen', 'green']  # 0.8, 0.9, 1.0

colors_85= ['mediumpurple', 'mediumorchid', 'plum', 'pink']  # 0.8, 0.9, 1.0
colors_86= ['gray', 'lightgray']
colors_88= ['darkred', 'firebrick', 'indianred', 'lightcoral']
colors_90= ['navy', 'tab:blue', 'skyblue', 'lightskyblue']  # 0.8, 0.9, 1.0
colors_95= ['green', 'limegreen', 'lawngreen', 'lightgreen']  # 0.8, 0.9, 1.0

color_tabl = [colors_85, colors_90, colors_95]


ylabels = {'cos_gwtheta': r'cos $\theta_\mathrm{gw}$',
           'gwphi': '$\phi_\mathrm{gw}$',
           'log10_h': 'log$_{10}h$',
           'log10_fgw': 'log$_{10}f_\mathrm{gw}$',
           'log10_Mc': 'log$_{10}M_\mathrm{c}$',
           'cos_inc': r'cos $i$',
           'phase0': r'$\Phi_0$',
           'psi': '$\Psi$'
           }


basepath ='./out/CGWparameter_recovery/ipta/'
folders = ['low', 'mid', 'high']

parameter = 'cos_gwtheta'


### 
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
plt.subplots_adjust(wspace=0.35)
for i, lmc in enumerate([8.5, 9.0, 9.5]):
    ivs = injected_values(22.3e-9, lmc)
    
    legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_tabl[i]]
    
    for j in range(3):
        data = np.genfromtxt(folders[j] + '/parameters_noPT_lmc{}_pd1.0.txt'.format(lmc), names=True)
        mask = np.argwhere((data['log10_fgw'] < -7.6) & (data['log10_fgw'] > -7.8))
        
        axs[i].errorbar(1+0.5*j*np.ones(len(data[parameter])), data[parameter], ls='', marker='x', color='k')        #yerr=[data[parameter]-data[parameter+'_ll'], data[parameter+'_ul']-data[parameter]]
        violin = axs[i].violinplot(data[parameter].transpose(), [1+0.5*j], showextrema=False, widths=0.2)
        for vl in violin['bodies']:
            vl.set_facecolor(color_tabl[i][j])
            vl.set_edgecolor(color_tabl[i][j])
            vl.set_alpha(0.5)
        
        
    axs[i].axhline(ivs[parameter], ls=':', color='gray')
    
    
    axs[i].set_title('log$_{10}M_\mathrm{c}$ = ' + '{}'.format(lmc))
    axs[i].set_xlabel('pulsar distance / kpc')
    axs[i].set_ylabel(ylabels[parameter])

plt.suptitle('Even Batches\n', fontsize=16, y=1.003)
plt.savefig('noPT_batches_{}.png'.format(parameter), bbox_inches='tight')
plt.show()





### dropout style ########################
folders = ['over1.0', 'over1.5', 'high']


### 
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
plt.subplots_adjust(wspace=0.35)
for i, lmc in enumerate([8.5, 9.0, 9.5]):
    ivs = injected_values(22.3e-9, lmc)
    
    legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_tabl[i]]
    
    for j in range(1,3):
        data = np.genfromtxt(folders[j] + '/parameters_noPT_lmc{}_pd1.0.txt'.format(lmc), names=True)
        mask = np.argwhere((data['log10_fgw'] < -7.6) & (data['log10_fgw'] > -7.8))
        
        axs[i].errorbar(1+0.5*j*np.ones(len(data[parameter])), data[parameter], ls='', marker='x', color='k')        #yerr=[data[parameter]-data[parameter+'_ll'], data[parameter+'_ul']-data[parameter]]
        violin = axs[i].violinplot(data[parameter].transpose(), [1+0.5*j], showextrema=False, widths=0.2)
        for vl in violin['bodies']:
            vl.set_facecolor(color_tabl[i][j])
            vl.set_edgecolor(color_tabl[i][j])
            vl.set_alpha(0.5)
        
        
    axs[i].axhline(ivs[parameter], ls=':', color='gray')
    
    
    axs[i].set_title('log$_{10}M_\mathrm{c}$ = ' + '{}'.format(lmc))
    axs[i].set_xlabel('pulsar distance / kpc')
    axs[i].set_ylabel(ylabels[parameter])

plt.suptitle('Dropout\n', fontsize=16, y=1.003)
plt.savefig('noPT_dropout_{}.png'.format(parameter), bbox_inches='tight')
plt.show()