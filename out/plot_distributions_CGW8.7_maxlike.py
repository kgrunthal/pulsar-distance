# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 14:41:08 2024

@author: kgrun
"""




import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches



plt.rcParams['font.family'] = "serif"
#plt.rcParams['font.sans-serif'] = "Times"

plt.rcParams['text.usetex']= False
plt.rcParams['xtick.labelsize'] = 10.0
plt.rcParams['ytick.labelsize'] = 10.0
plt.rcParams['axes.labelsize'] = 20.0



def plotconfig_new(fig, axs, fgw):
    for i in range(len(axs)):
        axs[i].axvline(fgw*1e9, ls=':', color='dimgray')
        
        axs[i].set_ylabel('$P(f)$', fontsize=14)
    
    axs[-1].set_xlabel('$f_\mathrm{gw}$ / nHz' , fontsize=14)
    axs[-1].set_ylabel('$P(f)$', fontsize=14)
    
    return None

fgw = 22.3e-9
bins = 20
test_frequencies = np.linspace(1/(3625*86400), 20/(3625*86400), bins)


realisations = 100

head_dir = './'


colors_85= ['mediumpurple', 'mediumorchid', 'plum']  # 0.8, 0.9, 1.0
colors_90= ['navy', 'tab:blue', 'skyblue']  # 0.8, 0.9, 1.0
colors_95= ['green', 'limegreen' , 'lawngreen']  # 0.8, 0.9, 1.0

color_tabl = [colors_85, colors_90, colors_95]

markers = ['D', 'o', '*', 'h']

legends = ['$d_\mathrm{p} = 1.0$kpc', '$d_\mathrm{p} = 4.0$kpc']


#### ISOTROPIC ################################################################
folders = ['WN_CGW_pd', 'isotropic_scaleto8.7', 'isotropic_5h']
lmcs = [8.7, 9.0, 8.7]
    

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(6, 14), sharex = True)
plt.subplots_adjust(wspace=0.02, hspace=0.1)

for ll, lmc in enumerate(lmcs):
    legend_patch= []
    
    for dd, pd in enumerate([1.0, 4.0]):
        OSarray = np.zeros((realisations, bins))
        SNarray = np.zeros((realisations, bins))
        for i in range(realisations):
            print('... on lmc {}, pd {}, realisation {}'.format(lmc, pd, i), end='\r')
            run_dir = folders[ll] + '/run_{}/'.format(i+1)
            
            test = np.loadtxt(run_dir + 'OS_spectrum_CGW{}_pd{}.txt'.format(lmc, pd)).transpose()
            OSarray[i] = test[1]
            SNarray[i] = test[1]/test[2]
            #axs[ll].errorbar(test_frequencies[1:9]*1e9+(dd*0.5-0.5), test[1][1:9]/test[2][1:9],
            #                 ls='', marker='.',color=color_tabl[ll][dd])
            axs[ll].errorbar(test_frequencies[1:9]*1e9+(dd*0.5-0.5), test[1][1:9],
                             ls='', marker='.',color=color_tabl[ll][dd])
        
        legend_patch.append(plt.Line2D((0,0), (0,0), ls='none', marker='.', color=color_tabl[ll][dd]))
    
    axs[ll].legend(legend_patch, legends, 
                   title = 'log$_{10}M_\mathrm{c}$ = ' + '{},  '.format(lmc) + '$h_0 = h_{0,(%d)}$'%(ll+1),
                   loc = 'upper left',
                   labelspacing=0.3, title_fontsize = 'large')

plotconfig_new(fig, axs, 22.3e-9)
plt.savefig(head_dir + 'isotropic-PFOS_CGW8.7_variations.png', dpi=400, bbox_inches='tight')
plt.show()
###############################################################################



#### GALACTIC #################################################################
folders = ['20PSR-galactic', 'galactic20_scaleto8.7', 'galactic20_5h']
lmcs = [8.7, 9.0, 8.7]
    

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(6, 14), sharex = True)
plt.subplots_adjust(wspace=0.02, hspace=0.1)

for ll, lmc in enumerate(lmcs):
    legend_patch= []
    
    for dd, pd in enumerate([1.0, 4.0]):
        OSarray = np.zeros((realisations, bins))
        SNarray = np.zeros((realisations, bins))
        for i in range(realisations):
            print('... on lmc {}, pd {}, realisation {}'.format(lmc, pd, i), end='\r')
            run_dir = folders[ll] + '/run_{}/'.format(i+1)
            
            test = np.loadtxt(run_dir + 'OS_spectrum_CGW{}_pd{}.txt'.format(lmc, pd)).transpose()
            OSarray[i] = test[1]
            axs[ll].errorbar(test_frequencies[1:9]*1e9+(dd*0.5-0.5), test[1][1:9],
                             ls='', marker='.',color=color_tabl[ll][dd])
           
        print()
        print(lmc, pd, np.argmax(OSarray[:,6]), np.max(OSarray[:,6])) 
        
        legend_patch.append(plt.Line2D((0,0), (0,0), ls='none', marker='.', color=color_tabl[ll][dd]))
    
    axs[ll].legend(legend_patch, legends, 
                   title = 'log$_{10}M_\mathrm{c}$ = ' + '{},  '.format(lmc) + '$h_0 = h_{0,(%d)}$'%(ll+1), 
                   loc = 'upper left',
                   labelspacing=0.3, title_fontsize = 'large')

plotconfig_new(fig, axs, 22.3e-9)
plt.savefig(head_dir + 'galactic20-PFOS_CGW8.7_variations.png', dpi=400, bbox_inches='tight')
plt.show()
###############################################################################