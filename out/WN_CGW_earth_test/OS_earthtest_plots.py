# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 11:51:08 2024

@author: kgrun
"""


import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


plt.rcParams['font.family'] = "serif"
plt.rcParams['font.sans-serif'] = "Times"

plt.rcParams['text.usetex']= False
plt.rcParams['xtick.labelsize'] = 11.0
plt.rcParams['ytick.labelsize'] = 11.0
plt.rcParams['axes.labelsize'] = 14.0



def plotaviolin(axs, data, OS_data, frequencies, facecolor, edgecolor, label=''):
    SN = []
    OS = []
    dOS = []
    
    OS_err = []
    OS_max, OS_range = [], []
    
    for i, f in enumerate(frequencies):
        OS.append(data[int(3*i)])
        dOS.append(data[int(3*i+1)])
        SN.append(data[int(3*i+2)])
        #SN.append(data[int(3*i)]/np.std(data[int(3*i)]))
        
        
        hist_OS, binedge_OS = np.histogram(data[int(3*i)], bins='auto')
        loc_OS = (binedge_OS[np.argmax(hist_OS)]+binedge_OS[np.argmax(hist_OS)+1])/2
        OS_max.append(loc_OS)
        OS_range.append(np.std(data[int(3*i)]))
        OS_err.append(np.mean(data[int(3*i+1)]))
        
    OS = np.array(OS)
    dOS = np.array(dOS)
    #SN = np.array(SN)
    SN = OS/dOS
    
    #axs[0].errorbar(np.log10(frequencies), OS_max, yerr=OS_range, ls='', fmt='rx', capsize=2)
    #axs[0].errorbar(np.log10(frequencies)+0.02  , OS_max, yerr=OS_err, ls='', fmt='bx', capsize=2)
    #print(OS_range[6],OS_err[6])
    #print(OS_range[0],OS_err[0])
    
    
    # left plot: OS values
    violin_OS = axs[0].violinplot(OS.transpose(), np.log10(frequencies),
                            showextrema=False, widths=0.05)
    for vl in violin_OS['bodies']:
        vl.set_facecolor(facecolor)
        vl.set_edgecolor(edgecolor)
        vl.set_alpha(0.5)
    
    # middle plot: dOS values
    violin_SN = axs[1].violinplot(dOS.transpose(), np.log10(frequencies),
                            showextrema=False, widths=0.05)
    for vl in violin_SN['bodies']:
        vl.set_facecolor(facecolor)
        vl.set_edgecolor(edgecolor)
        vl.set_alpha(0.5)
        
    # right plot: SN values
    violin_SN = axs[2].violinplot(SN.transpose(), np.log10(frequencies),
                            showextrema=False, widths=0.05)
    for vl in violin_SN['bodies']:
        vl.set_facecolor(facecolor)
        vl.set_edgecolor(edgecolor)
        vl.set_alpha(0.5)
    #plt.errorbar(np.log10(frequencies), OS_data[1]/OS_data[2],
    #             marker='o', markersize=4, ls='', color=edgecolor)

    legend = [mpatches.Patch(facecolor=facecolor, edgecolor=edgecolor, alpha=0.5), label]
    
    return legend



def plotconfig(fig, axs, fgw):
    #axs[1].yaxis.set_label_position("right")
    #axs[1].yaxis.tick_right()
    axs[2].set_ylim(-2)
    
    axs[0].axvline(np.log10(fgw), ls=':', color='darkgray')
    #axs[0].axvline(np.log10(0.8*fgw), ls=':', color='darkgray')
    axs[1].axvline(np.log10(fgw), ls=':', color='darkgray')
    #axs[1].axvline(np.log10(0.8*fgw), ls=':', color='darkgray')
    
    axs[0].set_xlabel('log$_{10}$($f$ / Hz)' , fontsize=14)
    axs[1].set_xlabel('log$_{10}$($f$ / Hz)' , fontsize=14)
    axs[2].set_xlabel('log$_{10}$($f$ / Hz)' , fontsize=14)
    
    #axs[0].set_ylabel('$\hat{A}_\mathrm{GW}^2$', fontsize=14)
    axs[0].set_ylabel('$P(f)$', fontsize=14)
    axs[1].set_ylabel('$\Delta_P(f)$', fontsize=14)
    axs[2].set_ylabel('S/N', fontsize=14)
    
    axs[0].set_title('PSD estimate distributions')
    axs[1].set_title('PSD $\sigma_0$ distributions')
    axs[2].set_title('S/N distributions')
    
    axs[0].set_yscale('log')
    axs[1].set_yscale('log')
    return None


    

def set_up_global_options():
    parser = argparse.ArgumentParser(description='Plot marginalised OS')
    parser.add_argument('--file', type=str, default=None, help='Path to the parfiles.')
    return parser.parse_args()



def main():
    args = set_up_global_options()
    
    fgw = 22.3e-9
    bins = 20
    test_frequencies = np.linspace(1/(3625*86400), 20/(3625*86400), bins)
    
    head_dir = './WN_CGW_earth_test/run_1/'
    ###########################################################################
    ###########################################################################
    outdir_maxLH = head_dir 
    
    OS_10 = np.loadtxt(outdir_maxLH + 'OS_spectrum_CGW10.0_pd1.0.txt')
    OS_95 = np.loadtxt(outdir_maxLH + 'OS_spectrum_CGW9.5_pd1.0.txt')
    OS_90 = np.loadtxt(outdir_maxLH + 'OS_spectrum_CGW9.0_pd1.0.txt')
    OS_85 = np.loadtxt(outdir_maxLH + 'OS_spectrum_CGW8.5_pd1.0.txt')
    
    ###########################################################################
    ###########################################################################
    outdir_nm = head_dir 
    
    data10 = np.loadtxt(outdir_nm + 'OS_spectrum_CGW10.0_pd1.0_NM.txt')
    data95 = np.loadtxt(outdir_nm + 'OS_spectrum_CGW9.5_pd1.0_NM.txt')
    data90 = np.loadtxt(outdir_nm + 'OS_spectrum_CGW9.0_pd1.0_NM.txt')
    data85 = np.loadtxt(outdir_nm + 'OS_spectrum_CGW8.5_pd1.0_NM.txt')
    
    ###########################################################################
    ###########################################################################

    
    
    outdir_figure = head_dir
    
    #-- logMc = 9.5 -----------------------------------------------------------
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
    plt.subplots_adjust(wspace=0.26)
    
    l95 = plotaviolin(axs, data95, OS_95, test_frequencies,
                      'lawngreen','limegreen',
                      label='log$M_\mathrm{c} = 9.5$')
    l90 = plotaviolin(axs, data90, OS_90, test_frequencies,
                      'skyblue','dodgerblue',
                      label='log$M_\mathrm{c} = 9.0$')
    l85 = plotaviolin(axs, data85, OS_85, test_frequencies,
                      'mediumpurple','rebeccapurple',
                      label='log$M_\mathrm{c} = 9.0$')
    
    plotconfig(fig, axs, fgw)
    plt.legend([l95[0], l90[0], l85[0]],
               [l95[1], l90[1], l85[1]],
               bbox_to_anchor = (1.6,1.))
    
    plt.savefig(outdir_figure + 'earthterm_OS_comparison.png', dpi=400, bbox_inches='tight')
    
    plt.show()
    
    
        
if __name__ == '__main__':
    main()