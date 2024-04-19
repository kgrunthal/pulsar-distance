# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 10:51:13 2024

@author: kgrun
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches



def plotaviolin(data, frequencies, facecolor, edgecolor, label=''):
    SN = []
    OS = []
    for i, f in enumerate(frequencies):
        SN.append(data[int(3*i+2)])
        OS.append(data[int(3*i)])
    SN = np.array(SN)  
    print(SN[1])
    OS = np.array(OS)
    violin = plt.violinplot(OS.transpose(), np.log10(frequencies),
                            showextrema=False, widths=0.05)
    
    #plt.errorbar(np.log10(frequencies), OS_data[1]/OS_data[2],
    #             marker='o', markersize=4, ls='', color=edgecolor)
    for vl in violin['bodies']:
        vl.set_facecolor(facecolor)
        vl.set_edgecolor(edgecolor)
        vl.set_alpha(0.5)
    legend = [mpatches.Patch(facecolor=facecolor, edgecolor=edgecolor, alpha=0.5), label]
    
    return violin, legend



OS_maxLH = np.loadtxt('./GWBsinglebin.txt').transpose()
OS_noisemargialised = np.loadtxt('./GWBsinglebin_NM.txt')

Tspan = 10*365.25*24*3600
freqs = (np.arange(30) + 1) / Tspan

#violin, legend = plotaviolin(OS_noisemargialised, freqs, 'lightblue', 'blue')
#plt.show()

fig,ax = plt.subplots()
ax_1 = ax.twinx()
ax.errorbar(np.log10(freqs), OS_maxLH[1], yerr=OS_maxLH[2], ls='', capsize=2, fmt='gx', label='OS')
ax_1.errorbar(np.log10(freqs), OS_maxLH[1]/OS_maxLH[2], ls='', fmt='bo', label='SN')
ax.set_xlabel('log10(frequency)')
ax.set_ylabel('OS value', color='green')
ax_1.set_ylabel('SN', color='blue')

fig.legend(bbox_to_anchor=(0.23,0.87))
plt.savefig('EPTA_OS_SN.png', bbox_inches='tight', dpi=400)
plt.show()
