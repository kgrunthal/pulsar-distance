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


colors_85= ['mediumpurple', 'mediumorchid', 'plum', 'pink']  # 0.8, 0.9, 1.0
colors_86= ['gray', 'lightgray']
colors_88= ['darkred', 'firebrick', 'indianred', 'lightcoral']
colors_90= ['navy', 'tab:blue', 'skyblue', 'lightskyblue']  # 0.8, 0.9, 1.0
colors_95= ['darkgreen', 'seagreen', 'limegreen', 'lightgreen']  # 0.8, 0.9, 1.0

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


basepath ='./'


parameter = 'cos_gwtheta'


#### BOXPLOTS ####


### BATCHES ###
subsets = ['low', 'mid', 'high']
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10, 4), sharey=True)
plt.subplots_adjust(wspace=0.05)
for i, lmc in enumerate([8.5, 9.0, 9.5]):
    ivs = injected_values(22.3e-9, lmc)
    
    legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_tabl[i]]
    
    for j, sub in enumerate(subsets):
        data = np.genfromtxt('parameters_noPT_lmc{}_{}.txt'.format(lmc, sub), names=True)
        
        boxprops = dict(facecolor=color_tabl[i][j],  edgecolor=color_tabl[i][j], linewidth=0, alpha=0.5)
        flierprops = dict(marker='o', markerfacecolor=color_tabl[i][j], ms=2.5, markeredgecolor='none')
        medianprops = dict(linewidth=1, color=color_tabl[i][j])
        whiskerprops = dict(color=color_tabl[i][j])
        capprops = dict(color=color_tabl[i][j])
        
        axs[i].boxplot(data[parameter].transpose(), positions=[1+0.5*j],
                          patch_artist=True,
                          boxprops=boxprops, flierprops=flierprops, medianprops=medianprops,
                          whiskerprops= whiskerprops, capprops=capprops)
        
        boxprops_1 = dict(color=color_tabl[i][j], linewidth=1)
        medianprops_1 = dict(linewidth=0)
        whiskerprops_1 = dict(color=color_tabl[i][j], linewidth=0)

        axs[i].boxplot(data[parameter].transpose(), positions=[1+0.5*j], patch_artist=False,
                          showmeans=False, showfliers=False, showcaps=False,
                          boxprops=boxprops_1, medianprops=medianprops_1,
                          whiskerprops= whiskerprops_1)
        
    axs[i].axhline(ivs[parameter], ls=':', color='gray')
    
    axs[i].set_xticks([1,1.5, 2])
    axs[i].set_xticklabels(['near', 'mid', 'far'])
    axs[i].set_title('log$_{10}M_\mathrm{c}$ = ' + '{}'.format(lmc))
    
    axs[i].set_xlim(0.8, 2.2)
axs[0].set_ylabel(ylabels[parameter])
    
        

#plt.suptitle('Even Batches\n', fontsize=16, y=1.003)
plt.savefig('noPT_batches_{}.png'.format(parameter), bbox_inches='tight', dpi=400)
plt.show()





### REMOVAL ########################
subsets = ['over1.0', 'over1.5', 'over2.0']

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10, 4), sharey=True)
plt.subplots_adjust(wspace=0.05)
for i, lmc in enumerate([8.5, 9.0, 9.5]):
    ivs = injected_values(22.3e-9, lmc)
    
    legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_tabl[i]]
    
    for j, sub in enumerate(subsets):
        data = np.genfromtxt('parameters_noPT_lmc{}_{}.txt'.format(lmc, sub), names=True)
        
        
        boxprops = dict(facecolor=color_tabl[i][j],  edgecolor=color_tabl[i][j], linewidth=0, alpha=0.5)
        flierprops = dict(marker='o', markerfacecolor=color_tabl[i][j], ms=2.5, markeredgecolor='none')
        medianprops = dict(linewidth=1, color=color_tabl[i][j])
        whiskerprops = dict(color=color_tabl[i][j])
        capprops = dict(color=color_tabl[i][j])
        
        axs[i].boxplot(data[parameter].transpose(), positions=[1+0.5*j],
                          patch_artist=True,
                          boxprops=boxprops, flierprops=flierprops, medianprops=medianprops,
                          whiskerprops= whiskerprops, capprops=capprops)
        
        boxprops_1 = dict(color=color_tabl[i][j], linewidth=1)
        medianprops_1 = dict(linewidth=0)
        whiskerprops_1 = dict(color=color_tabl[i][j], linewidth=0)

        axs[i].boxplot(data[parameter].transpose(), positions=[1+0.5*j], patch_artist=False,
                          showmeans=False, showfliers=False, showcaps=False,
                          boxprops=boxprops_1, medianprops=medianprops_1,
                          whiskerprops= whiskerprops_1)
        
    axs[i].axhline(ivs[parameter], ls=':', color='gray')
    
    axs[i].set_xticks([1,1.5, 2])
    axs[i].set_xticklabels(['> 1 kpc', '> 1.5 kpc', '> 2 kpc'])
    axs[i].set_title('log$_{10}M_\mathrm{c}$ = ' + '{}'.format(lmc))
    
    axs[i].set_xlim(0.8, 2.2)

axs[0].set_ylabel(ylabels[parameter])
    
    
#plt.suptitle('Dropout\n', fontsize=16, y=1.003)
plt.savefig('noPT_removal_{}.png'.format(parameter), bbox_inches='tight', dpi=400)
plt.show()


'''
### BATCHES TEST ########################
folders = ['mid_as_low', 'mid', 'mid_as_high']



fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10, 4), sharey=True)
plt.subplots_adjust(wspace=0.05)
for i, lmc in enumerate([8.5, 9.0, 9.5]):
    ivs = injected_values(22.3e-9, lmc)
    
    legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_tabl[i]]
    
    for j in range(0,3):
        data = np.genfromtxt(folders[j] + '/parameters_noPT_lmc{}_pd1.0.txt'.format(lmc), names=True)
        
        boxprops = dict(facecolor=color_tabl[i][j],  edgecolor=color_tabl[i][j], linewidth=0, alpha=0.5)
        flierprops = dict(marker='o', markerfacecolor=color_tabl[i][j], ms=2.5, markeredgecolor='none')
        medianprops = dict(linewidth=1, color=color_tabl[i][j])
        whiskerprops = dict(color=color_tabl[i][j])
        capprops = dict(color=color_tabl[i][j])
        
        axs[i].boxplot(data[parameter].transpose(), positions=[1+0.5*j],
                          patch_artist=True,
                          boxprops=boxprops, flierprops=flierprops, medianprops=medianprops,
                          whiskerprops= whiskerprops, capprops=capprops)
        
        boxprops_1 = dict(color=color_tabl[i][j], linewidth=1)
        medianprops_1 = dict(linewidth=0)
        whiskerprops_1 = dict(color=color_tabl[i][j], linewidth=0)

        axs[i].boxplot(data[parameter].transpose(), positions=[1+0.5*j], patch_artist=False,
                          showmeans=False, showfliers=False, showcaps=False,
                          boxprops=boxprops_1, medianprops=medianprops_1,
                          whiskerprops= whiskerprops_1)
        
    axs[i].axhline(ivs[parameter], ls=':', color='gray')
    
    
    axs[i].set_xticks([1,1.5, 2])
    axs[i].set_xticklabels(['$-$1 kpc', '+ 0 kpc', '$+$ 1 kpc'])
    axs[i].set_title('log$_{10}M_\mathrm{c}$ = ' + '{}'.format(lmc))
    axs[i].set_xlabel('mid dataset')
    
    axs[i].set_xlim(0.8, 2.2)
axs[0].set_ylabel(ylabels[parameter])

    
    
#plt.suptitle('Dropout\n', fontsize=16, y=1.003)
plt.savefig('noPT_batchestest_{}.png'.format(parameter), bbox_inches='tight', dpi=400)
plt.show()

'''







