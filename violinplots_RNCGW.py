# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 12:30:17 2024

@author: kgrun
"""


import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches



def plotaviolin(data, OS_data, frequencies, facecolor, edgecolor, label=''):
    SN = []
    OS = []
    for i, f in enumerate(frequencies):
        SN.append(data[int(3*i+2)])
        OS.append(data[int(3*i)])
    SN = np.array(SN)  
    OS = np.array(OS)
    
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(6, 3))
    
    # left plot: OS values
    violin_OS = axs[0].violinplot(OS.transpose(), np.log10(frequencies),
                            showextrema=False, widths=0.05)
    for vl in violin_OS['bodies']:
        vl.set_facecolor(facecolor)
        vl.set_edgecolor(edgecolor)
        vl.set_alpha(0.5)
    
    # right plot: SN values
    violin_SN = axs[1].violinplot(SN.transpose(), np.log10(frequencies),
                            showextrema=False, widths=0.05)
    for vl in violin_SN['bodies']:
        vl.set_facecolor(facecolor)
        vl.set_edgecolor(edgecolor)
        vl.set_alpha(0.5)
    #plt.errorbar(np.log10(frequencies), OS_data[1]/OS_data[2],
    #             marker='o', markersize=4, ls='', color=edgecolor)
    
    legend = [mpatches.Patch(facecolor=facecolor, edgecolor=edgecolor, alpha=0.5), label]
    
    return legend


def set_up_global_options():
    parser = argparse.ArgumentParser(description='Plot marginalised OS')
    parser.add_argument('--file', type=str, default=None, help='Path to the parfiles.')
    return parser.parse_args()



def main():
    args = set_up_global_options()
    
    fgw = 22.3e-9
    bins = 20
    test_frequencies = np.linspace(1/(3625*86400), 20/(3625*86400), bins)
    
    
    OS_90_10 = np.loadtxt('./out/OS_spectrum_RNCGW9.0_zeta1.0.txt')
    data90_10 = np.loadtxt('./out/OS_spectrum_RNCGW9.0_zeta1.0_NM.txt')
    OS_90_10_WN = np.loadtxt('./out/maxLH/OS_spectrum_WNCGW9.0_zeta1.0.txt').transpose()
    data90_10_WN = np.loadtxt('./out/noisemarginalised/OS_spectrum_WNCGW9.0_zeta1.0_NM.txt').transpose()
    
    
    plt.figure(figsize=(8,6))
    
    violin90_10, l90_10 = plotaviolin(data90_10, OS_90_10, test_frequencies,
                                      'navy','midnightblue',
                                      label='log$M_\mathrm{c} = 9.0$, $\zeta = 1.0$')
    violin90_10_WN, l90_10_WN = plotaviolin(data90_10_WN, OS_90_10_WN, test_frequencies,
                                      'green','darkgreen',
                                      label='log$M_\mathrm{c} = 9.0$, $\zeta = 1.0$')
    plt.axvline(np.log10(fgw), ls=':', color='darkgray')
    
    plt.xlabel('log$_{10}$(frequency / Hz)' , fontsize=14)
    plt.ylabel('S/N', fontsize=14)
    plt.title('log Mc = 9.0')
    #plt.legend([l90_08[0], l90_09[0], l90_10[0]],
    #           [l90_08[1], l90_09[1], l90_10[1]])
    plt.legend([l90_10[0]], [l90_10[1]])
    plt.show()
    
    
if __name__ == '__main__':
    main()