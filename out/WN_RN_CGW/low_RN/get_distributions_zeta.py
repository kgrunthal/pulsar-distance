# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 13:14:25 2024

@author: kgrun
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches

def gen_color(cmap,n,reverse=False):
    '''Generates n distinct color from a given colormap.

    Args:
        cmap(str): The name of the colormap you want to use.
            Refer https://matplotlib.org/stable/tutorials/colors/colormaps.html to choose
            Suggestions:
            For Metallicity in Astrophysics: Use coolwarm, bwr, seismic in reverse
            For distinct objects: Use gnuplot, brg, jet,turbo.

        n(int): Number of colors you want from the cmap you entered.

        reverse(bool): False by default. Set it to True if you want the cmap result to be reversed.

    Returns: 
        colorlist(list): A list with hex values of colors.
    '''
    c_map = plt.cm.get_cmap(str(cmap)) # select the desired cmap
    arr=np.linspace(0,1,n) #create a list with numbers from 0 to 1 with n items
    colorlist=list()
    for c in arr:
        rgba=c_map(c) #select the rgba value of the cmap at point c which is a number between 0 to 1
        clr=colors.rgb2hex(rgba) #convert to hex
        colorlist.append(str(clr)) # create a list of these colors
    
    if reverse==True:
        colorlist.reverse()
    return colorlist


def plotconfig(fig, axs, fgw, style='log'):
    if style == 'log':
        axs[0].axvline(np.log10(fgw), ls=':', color='dimgray')
        axs[0].axvline(np.log10(0.9*fgw), ls=':', color='darkgray')
        axs[0].axvline(np.log10(0.8*fgw), ls=':', color='darkgray')
        
        axs[1].axvline(np.log10(fgw), ls=':', color='dimgray')
        axs[1].axvline(np.log10(0.9*fgw), ls=':', color='darkgray')
        axs[1].axvline(np.log10(0.8*fgw), ls=':', color='darkgray')
    
        axs[0].set_xlabel('log$_{10}$($f$ / Hz)' , fontsize=14)
        axs[1].set_xlabel('log$_{10}$($f$ / Hz)' , fontsize=14)
        axs[0].set_ylabel('$\hat{A}_\mathrm{GW}^2$', fontsize=14)
        axs[1].set_ylabel('S/N', fontsize=14)
    
    elif style == 'even':
        axs[0].axvline(fgw, ls=':', color='dimgray')
        axs[0].axvline(0.9*fgw, ls=':', color='darkgray')
        axs[0].axvline(0.8*fgw, ls=':', color='darkgray')
        
        bin_pos = np.arange(0, 21, 5)
        axs[0].set_xticks(bin_pos)
        axs[0].set_xlim(0,11.5)
        axs[1].set_xticks(bin_pos)
        axs[1].set_xlim(0,11.5)
        axs[2].set_xticks(bin_pos)
        axs[2].set_xlim(0,21)
        #axs[1].set_xlim(right=3.5e-8)
        
        axs[0].set_xlabel('bin number' , fontsize=14)
        axs[1].set_xlabel('bin number' , fontsize=14)
        axs[2].set_xlabel('bin number' , fontsize=14)
        axs[0].set_ylabel('$P(f)$', fontsize=14)
        axs[1].set_ylabel('$\Delta_P(f)$', fontsize=14)
        axs[2].set_ylabel('S/N', fontsize=14)
    
    
    return None


def distribution(data, frequencies):
    SN_max = []
    OS_max = []
    
    for i, f in enumerate(frequencies):
        OS = data[int(3*i)]*f
        SN = data[int(3*i+2)]
        
        
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
bin_numbers = np.arange(1,21,step=1)


head_dir = './out/WN_RN_CGW/'

realisations = 10

color_08 = gen_color(cmap='Greens', n=20, reverse=True)
color_09 = gen_color(cmap='Blues', n=20, reverse=True)
color_10 = gen_color(cmap='Purples', n=20, reverse=True)


cmaps = [color_08, color_09, color_10]

colors_85= ['plum', 'mediumorchid', 'mediumpurple']  # 0.8, 0.9, 1.0
colors_90= ['skyblue', 'tab:blue', 'navy']  # 0.8, 0.9, 1.0
colors_95= ['lawngreen', 'limegreen', 'green']  # 0.8, 0.9, 1.0
color_tabl = [colors_85, colors_90, colors_95]

markers = ['D', 'o', '*', 'h']

legends = ['$\zeta = 0.8$', '$\zeta = 0.9$', '$\zeta = 1.0$']

'''
for ll, lmc in enumerate([8.5, 9.0, 9.5]):
    
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(20, 5))
    plt.subplots_adjust(wspace=0.05)
    axs[1].yaxis.set_label_position("right")
    axs[1].yaxis.tick_right() 
    plotconfig(fig, axs, fgw)
    
    legend_patch= []
    for zz, zeta in enumerate([0.8, 0.9, 1.0]):
        OSarray = np.zeros((realisations, bins))
        SNarray = np.zeros((realisations, bins))
        for i in range(realisations):
            run_dir = 'run_{}/noisemarginalised/'.format(i+1)
            data = np.loadtxt(head_dir + run_dir + 'OS_spectrum_CGW{}_zeta{}_NM.txt'.format(lmc, zeta))
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
        
        axs[0].errorbar(np.log10(test_frequencies)+zz*0.01, OS_avg, yerr = OS_range,
                       ls='', markersize=0, capsize=3, linewidth=3, elinewidth=3,
                       color=color_tabl[ll][zz])
        axs[1].errorbar(np.log10(test_frequencies)+zz*0.01, SN_avg, yerr = SN_range,
                       ls='', markersize=0, capsize=3, linewidth=3, elinewidth=3,
                       color=color_tabl[ll][zz])
        legend_patch.append(plt.Line2D((0,0), (0,0), ls='none', marker=markers[3], color=color_tabl[ll][zz]))
    
    
    
    plt.suptitle('log$M_\mathrm{c} =$' + '{}'.format(lmc), fontsize=20)
    plt.legend(legend_patch, legends, loc = 'upper left')
    #plt.savefig(head_dir + 'realisations_{}.png'.format(lmc), dpi=400, bbox_inches='tight')
    plt.show()
'''        
        

for ll, lmc in enumerate([8.5, 9.0, 9.5]):
    
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(20, 5))
    plotconfig(fig, axs, fgw, style='even')
    
    legend_patch= []
    off=0.2
    for zz, zeta in enumerate([0.8, 0.9, 1.0]):
        OSarray = np.zeros((realisations, bins))
        SNarray = np.zeros((realisations, bins))
        for i in range(realisations):
            run_dir = 'run_{}/noisemarginalised/'.format(i+1)
            data = np.loadtxt(head_dir + run_dir + 'OS_spectrum_RNCGW{}_zeta{}_NM.txt'.format(lmc, zeta))
            OS, SN = distribution(data, test_frequencies)
            OSarray[i], SNarray[i] = OS, SN
            
            axs[0].errorbar(bin_numbers+zz*off, OS, ls='',
                         color=color_tabl[ll][zz], marker=markers[zz], ms=2)
            
            axs[1].errorbar(bin_numbers+zz*off, SN, ls='',
                         color=color_tabl[ll][zz], marker=markers[zz], ms=2)
        
        # produce ranges for plotting
        OS_avg = (np.max(OSarray, axis=0)+ np.min(OSarray, axis=0))/2.
        OS_range = np.max(OSarray, axis=0) - OS_avg
        SN_avg = (np.max(SNarray, axis=0)+np.min(SNarray, axis=0))/2
        SN_range = np.max(SNarray, axis=0) - SN_avg
        
        axs[0].errorbar(bin_numbers+zz*off, OS_avg, yerr = OS_range,
                       ls='', markersize=0, capsize=0, linewidth=3, elinewidth=5,
                       color=color_tabl[ll][zz], alpha=0.5)
        
        
        axs[1].errorbar(bin_numbers + zz*off, SN_avg, yerr = SN_range,
                       ls='', markersize=0, capsize=0, linewidth=3, elinewidth=5,
                       color=color_tabl[ll][zz], alpha=0.5)
        axs[2].errorbar(bin_numbers, OS_avg/OS_range,
                        ls='', linewidth=3, elinewidth=5, marker=markers[zz],
                        color=color_tabl[ll][zz])
        legend_patch.append(plt.Line2D((0,0), (0,0), ls='none', marker=markers[zz], color=color_tabl[ll][zz]))
    
    
    
    plt.suptitle('log$M_\mathrm{c} =$' + '{}'.format(lmc), fontsize=20)
    plt.legend(legend_patch, legends, loc = 'upper left')
    #plt.savefig(head_dir + 'realisations_{}.png'.format(lmc), dpi=400, bbox_inches='tight')
    plt.show()
        


