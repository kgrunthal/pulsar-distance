# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 13:14:25 2024

@author: kgrun
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches


def plotconfig(fig, axs, fgw):
    axs[0].axvline(np.log10(fgw), ls=':', color='dimgray')
    axs[0].axvline(np.log10(0.9*fgw), ls=':', color='darkgray')
    axs[0].axvline(np.log10(0.8*fgw), ls=':', color='darkgray')
    axs[1].axvline(np.log10(fgw), ls=':', color='dimgray')
    axs[1].axvline(np.log10(0.9*fgw), ls=':', color='darkgray')
    axs[1].axvline(np.log10(0.8*fgw), ls=':', color='darkgray')
    
    axs[0].set_xlabel('log$_{10}$($f$ / Hz)' , fontsize=14)
    axs[1].set_xlabel('log$_{10}$($f$ / Hz)' , fontsize=14)
    axs[0].set_ylabel('$P(f)$', fontsize=14)
    axs[1].set_ylabel(r'$\overline{P} / (P_\mathrm{max}-P_\mathrm{min})$', fontsize=14)
    
    return None


def distribution(data, frequencies):
    SN_max = []
    OS_max = []
    
    for i, f in enumerate(frequencies):
        SN = data[int(3*i+2)]
        OS = data[int(3*i)]
        
        hist_SN, binedge_SN = np.histogram(SN, bins='auto')
        loc_SN = (binedge_SN[np.argmax(hist_SN)]+binedge_SN[np.argmax(hist_SN)+1])/2
        SN_max.append(loc_SN)
        
        hist_OS, binedge_OS = np.histogram(OS, bins='auto')
        loc_OS = (binedge_OS[np.argmax(hist_OS)]+binedge_OS[np.argmax(hist_OS)+1])/2
        OS_max.append(loc_OS)
        
    SN_max = np.array(SN_max)  
    OS_max = np.array(OS_max)
    
    return OS_max, SN_max
    
    

fgw = 22.3e-9
bins = 20
test_frequencies = np.linspace(1/(3625*86400), 20/(3625*86400), bins)


realisations = 50



colors_85= ['mediumpurple', 'mediumorchid', 'plum']  # 0.8, 0.9, 1.0
colors_90= ['navy', 'tab:blue', 'skyblue']  # 0.8, 0.9, 1.0
colors_95= ['green', 'limegreen' , 'lawngreen']  # 0.8, 0.9, 1.0

color_tabl = [colors_85, colors_90, colors_95]

markers = ['D', 'o', '*', 'h']

legends = ['$\zeta = 1.0$', '$\zeta = 0.9$', '$\zeta = 0.8$']


'''
for ll, lmc in enumerate([8.5, 9.0, 9.5]):
    
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(20, 5))
    plt.subplots_adjust(wspace=0.05)
    axs[1].yaxis.set_label_position("right")
    axs[1].yaxis.tick_right() 
    plotconfig(fig, axs, fgw)
    
    legend_patch= []
    for dd, pd in enumerate([1.0, 1.5, 2.0]):
        OSarray = np.zeros((realisations, bins))
        SNarray = np.zeros((realisations, bins))
        for i in range(realisations):
            run_dir = 'run_{}/noisemarginalised/'.format(i+1)
            data = np.loadtxt(head_dir + run_dir + 'OS_spectrum_CGW{}_pd{}_NM.txt'.format(lmc, pd))
            OS, SN = distribution(data, test_frequencies)
            OSarray[i], SNarray[i] = OS, SN
            
            #axs[0].errorbar(np.log10(test_frequencies), OS, ls='',
            #             color=color_tabl[ll][zz], marker=markers[2])
            
            #axs[1].errorbar(np.log10(test_frequencies), SN, ls='',
            #             color=color_tabl[ll][zz], marker=markers[2])
        
        # produce ranges for plotting
        OS_avg = np.average(OSarray, axis=0)
        OS_range = np.max(OSarray, axis=0) - OS_avg
        SN_avg = np.average(SNarray, axis=0)
        SN_range = np.max(SNarray, axis=0) - SN_avg
        
        axs[0].errorbar(np.log10(test_frequencies)+dd*0.01, OS_avg, yerr = OS_range,
                       ls='', markersize=0, capsize=3, linewidth=3, elinewidth=3,
                       color=color_tabl[ll][dd])
        axs[1].errorbar(np.log10(test_frequencies)+dd*0.01, SN_avg, yerr = SN_range,
                       ls='', markersize=0, capsize=3, linewidth=3, elinewidth=3,
                       color=color_tabl[ll][dd])
        legend_patch.append(plt.Line2D((0,0), (0,0), ls='none', marker=markers[3], color=color_tabl[ll][dd]))
    
    
    
    plt.suptitle('log$M_\mathrm{c} =$' + '{}'.format(lmc), fontsize=20)
    plt.legend(legend_patch, legends, loc = 'upper left')
    plt.savefig(head_dir + 'realisations_{}.png'.format(lmc), dpi=400, bbox_inches='tight')
    plt.show()
'''

'''
for ll, lmc in enumerate([8.5, 9.0, 9.5]):
    
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
    plt.subplots_adjust(wspace=0.02)
    axs[1].yaxis.set_label_position("right")
    axs[1].yaxis.tick_right() 
    plotconfig(fig, axs, fgw)
    
    legend_patch= []
    for dd, zeta in enumerate([1.0, 0.9, 0.8]):
        OSarray = np.zeros((realisations, bins))
        SNarray = np.zeros((realisations, bins))
        for i in range(realisations):
            run_dir = 'run_{}/noisemarginalised/'.format(i+1)
            data = np.loadtxt(run_dir + 'OS_spectrum_CGW{}_zeta{}_NM.txt'.format(lmc, zeta))
            OS, SN = distribution(data, test_frequencies)
            OSarray[i], SNarray[i] = OS, SN
            
            #axs[0].errorbar(np.log10(test_frequencies)+dd*0.01, OS, ls='',
            #            color=color_tabl[ll][dd], marker=markers[2], markersize=3)
            
            #axs[1].errorbar(np.log10(test_frequencies)+dd*0.01, SN, ls='',
            #             color=color_tabl[ll][dd], marker=markers[2], markersize=3)
        
        # produce ranges for plotting
        OS_avg = np.average(OSarray, axis=0)
        OS_range = np.max(OSarray, axis=0) - OS_avg
        SN_avg = np.average(SNarray, axis=0)
        SN_range = np.max(SNarray, axis=0) - SN_avg
        
        
        # old verison with OS/dOS
        #axs[0].errorbar(np.log10(test_frequencies)+dd*0.01, OS_avg, yerr = OS_range,
        #               ls='', markersize=0, capsize=3, linewidth=3, elinewidth=3,
        #               color=color_tabl[ll][dd])
        #axs[1].errorbar(np.log10(test_frequencies)+dd*0.01, SN_avg, yerr = SN_range,
        #               ls='', markersize=0, capsize=3, linewidth=3, elinewidth=3,
        #               color=color_tabl[ll][dd])
            
        
        
        # new version with spread of distribution
        axs[0].errorbar(np.log10(test_frequencies)+dd*0.015, OS_avg, yerr = OS_range,
                       ls='', markersize=0, capsize=3, linewidth=3, elinewidth=3,
                       color=color_tabl[ll][dd])
        
        axs[1].errorbar(np.log10(test_frequencies), OS_avg/OS_range,
                        ls='', linewidth=3, elinewidth=3, marker=markers[dd],
                        color=color_tabl[ll][dd])
        
        legend_patch.append(plt.Line2D((0,0), (0,0), ls='none', marker=markers[dd], color=color_tabl[ll][dd]))
        
        
        # violin plots
        #axs[0].violinplot(OSarray, np.log10(test_frequencies),showextrema=False, widths=0.05)
        #axs[1].violin(np.log10(test_frequencies), SNarray.transpose())
    
    
    #plt.suptitle('log$M_\mathrm{c} =$' + '{}'.format(lmc), fontsize=20)
    plt.legend(legend_patch, legends, loc = 'upper left')
    plt.savefig('PFOS_realisations_{}_new.png'.format(lmc), dpi=400, bbox_inches='tight')
    plt.show()
'''       




#### ALL TOGETHER IN ONE ####
def plotconfig_new(fig, axs, fgw):
    for i in range(3):
        axs[i][0].axvline(fgw*1e9, ls=':', color='dimgray')
        axs[i][0].axvline(0.9*fgw*1e9, ls=':', color='darkgray')
        axs[i][0].axvline(0.8*fgw*1e9, ls=':', color='darkgray')
        axs[i][1].axvline(fgw*1e9, ls=':', color='dimgray')
        axs[i][1].axvline(0.9*fgw*1e9, ls=':', color='darkgray')
        axs[i][1].axvline(0.8*fgw*1e9, ls=':', color='darkgray')
        
        axs[i][0].set_ylabel('$P(f)$', fontsize=14)
        axs[i][1].set_ylabel(r'$\overline{P} / (P_\mathrm{max}-P_\mathrm{min})$', fontsize=14)
    
    axs[2][0].set_xlabel('$f_\mathrm{gw}$ / nHz' , fontsize=14)
    axs[2][1].set_xlabel('$f_\mathrm{gw}$ / nHz' , fontsize=14)
    axs[2][0].set_ylabel('$P(f)$', fontsize=14)
    axs[2][1].set_ylabel(r'$\overline{P} / (P_\mathrm{max}-P_\mathrm{min})$', fontsize=14)
    
    return None
        


fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(10, 13), sharex = True)
plt.subplots_adjust(wspace=0.02, hspace=0.1)
axs[0][1].yaxis.set_label_position("right")
axs[0][1].yaxis.tick_right() 

for ll, lmc in enumerate([8.5, 9.0, 9.5]):
    print('On lmc {}'.format(lmc))
    
    axs[ll][1].yaxis.set_label_position("right")
    axs[ll][1].yaxis.tick_right() 
    legend_patch= []
    
    for dd, zeta in enumerate([1.0, 0.9, 0.8]):
        print('... zeta {} \t'.format(zeta), end='\r')
        OSarray = np.zeros((realisations, bins))
        SNarray = np.zeros((realisations, bins))
        for i in range(realisations):
            run_dir = 'run_{}/noisemarginalised/'.format(i+1)
            data = np.loadtxt(run_dir + 'OS_spectrum_CGW{}_zeta{}_NM.txt'.format(lmc, zeta))
            OS, SN = distribution(data, test_frequencies)
            OSarray[i], SNarray[i] = OS, SN
        print('DONE')
        # produce ranges for plotting
        OS_avg = np.average(OSarray, axis=0)
        OS_range = np.max(OSarray, axis=0) - OS_avg
        SN_avg = np.average(SNarray, axis=0)
        SN_range = np.max(SNarray, axis=0) - SN_avg
        
       
        # new version with spread of distribution
        axs[ll][0].errorbar(test_frequencies*1e9+dd*0.5, OS_avg, yerr = OS_range,
                       ls='', markersize=0, capsize=3, linewidth=3, elinewidth=3,
                       color=color_tabl[ll][dd])
        
        axs[ll][1].errorbar(test_frequencies*1e9, OS_avg/OS_range,
                        ls='', linewidth=3, elinewidth=3, marker=markers[dd],
                        color=color_tabl[ll][dd], label='$\zeta = {}$'.format(zeta))
        
        legend_patch.append(plt.Line2D((0,0), (0,0), ls='none', marker=markers[dd], color=color_tabl[ll][dd]))
        
          
    
    
    axs[ll][1].legend(loc = 'upper right')

plotconfig_new(fig, axs, fgw)

plt.savefig('ring-PFOS_realisations.png', dpi=400, bbox_inches='tight')
plt.show()

