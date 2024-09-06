# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 09:45:31 2024

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
        OS.append(data[int(3*i)]*f)
        dOS.append(data[int(3*i+1)]*f)
        SN.append(data[int(3*i+2)])
        #SN.append(data[int(3*i)]/np.std(data[int(3*i)]))
        
        
        hist_OS, binedge_OS = np.histogram(data[int(3*i)], bins='auto')
        loc_OS = (binedge_OS[np.argmax(hist_OS)]+binedge_OS[np.argmax(hist_OS)+1])/2
        OS_max.append(loc_OS)
        OS_range.append(np.std(data[int(3*i)]))
        OS_err.append(np.mean(data[int(3*i+1)]))
        
       
    OS = np.array(OS)
    dOS = np.array(dOS)
    SN = np.array(SN)
    #SN = OS/dOS
    
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
    axs[2].set_ylim(-4)
    
    axs[0].axvline(np.log10(fgw), ls=':', color='darkgray')
    axs[0].axvline(np.log10(0.8*fgw), ls=':', color='darkgray')
    axs[1].axvline(np.log10(fgw), ls=':', color='darkgray')
    axs[1].axvline(np.log10(0.8*fgw), ls=':', color='darkgray')
    
    axs[0].set_xlabel('log$_{10}$($f$ / Hz)' , fontsize=14)
    axs[1].set_xlabel('log$_{10}$($f$ / Hz)' , fontsize=14)
    axs[2].set_xlabel('log$_{10}$($f$ / Hz)' , fontsize=14)
    
    #axs[0].set_ylabel('$\hat{A}_\mathrm{GW}^2$', fontsize=14)
    axs[0].set_ylabel('$P(f)$', fontsize=14)
    axs[1].set_ylabel('$\Delta_{P(f)}$', fontsize=14)
    axs[2].set_ylabel('S/N', fontsize=14)
    
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
    
    
    ###########################################################################
    ###########################################################################
    head_dir = './out/WN_RN_CGW_higher/'
    
    ###########################################################################
    ###########################################################################
    outdir_maxLH = head_dir + '/maxLH/'
    
    OS_95_08 = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW9.5_zeta0.8.txt')
    OS_95_09 = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW9.5_zeta0.9.txt')
    OS_95_10 = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW9.5_zeta1.0.txt')
    #OS_95_earth = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW9.5_earth.txt')
    #OS_95_pulsar = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW9.5_pulsar.txt')
    
    OS_90_08 = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW9.0_zeta0.8.txt')
    OS_90_09 = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW9.0_zeta0.9.txt')
    OS_90_10 = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW9.0_zeta1.0.txt')
    #OS_90_earth = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW9.0_earth.txt')
    #OS_90_pulsar = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW9.0_pulsar.txt')
    
    OS_85_08 = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW8.5_zeta0.8.txt')
    OS_85_09 = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW8.5_zeta0.9.txt')
    OS_85_10 = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW8.5_zeta1.0.txt')
    #OS_85_earth = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW8.5_earth.txt')
    #OS_85_pulsar = np.loadtxt(outdir_maxLH + 'OS_spectrum_RNCGW8.5_pulsar.txt')
    
    #OS_WNGWB = np.loadtxt(outdir_maxLH + 'OS_spectrum_WNGWB.txt').transpose()
    
    ###########################################################################
    ###########################################################################
    outdir_nm = head_dir + '/noisemarginalised/'
    
    data95_08 = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW9.5_zeta0.8_NM.txt')
    data95_09 = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW9.5_zeta0.9_NM.txt')
    data95_10 = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW9.5_zeta1.0_NM.txt')
    #data95_earth = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW9.5_earth_NM.txt')
    #data95_pulsar = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW9.5_pulsar_NM.txt')
    
    data90_08 = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW9.0_zeta0.8_NM.txt')
    data90_09 = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW9.0_zeta0.9_NM.txt')
    data90_10 = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW9.0_zeta1.0_NM.txt')
    #data90_earth = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW9.0_earth_NM.txt')
    #data90_pulsar = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW9.0_pulsar_NM.txt')
    
    data85_08 = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW8.5_zeta0.8_NM.txt')
    data85_09 = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW8.5_zeta0.9_NM.txt')
    data85_10 = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW8.5_zeta1.0_NM.txt')
    #data85_earth = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW8.5_earth_NM.txt')
    #data85_pulsar = np.loadtxt(outdir_nm + 'OS_spectrum_RNCGW8.5_pulsar_NM.txt')
    
    ###########################################################################
    ###########################################################################
    outdir_figure = head_dir
    
    #-- logMc = 9.5 -----------------------------------------------------------
    fig_95, axs_95 = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
    #plt.subplots_adjust(wspace=0.05)
    
    l95_10 = plotaviolin(axs_95, data95_10, OS_95_10, test_frequencies,
                         'green','darkgreen',
                         label='log$M_\mathrm{c} = 9.5$, $\zeta = 1.0$')
    l95_09 = plotaviolin(axs_95, data95_09, OS_95_09, test_frequencies,
                         'limegreen','forestgreen',
                         label='log$M_\mathrm{c} = 9.5$, $\zeta = 0.9$')
    l95_08 = plotaviolin(axs_95, data95_08, OS_95_08, test_frequencies,
                         'lawngreen','limegreen',
                         label='log$M_\mathrm{c} = 9.5$, $\zeta = 0.8$')
    plotconfig(fig_95, axs_95, fgw)
   
    plt.suptitle('log Mc = 9.5', fontsize=20)
    plt.legend([l95_08[0], l95_09[0], l95_10[0]],
               [l95_08[1], l95_09[1], l95_10[1]], 
               bbox_to_anchor = (1.5,1.))
    plt.savefig(outdir_figure + 'logMc9.5.png', dpi=400, bbox_inches='tight')
    plt.show()
    #--------------------------------------------------------------------------
    
    #-- logMc = 9.0 -----------------------------------------------------------
    fig_90, axs_90 = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
    #plt.subplots_adjust(wspace=0.05)
    
    l90_08 = plotaviolin(axs_90, data90_08, OS_90_08, test_frequencies,
                         'skyblue','dodgerblue',
                         label='log$M_\mathrm{c} = 9.0$, $\zeta = 0.8$')
    l90_09 = plotaviolin(axs_90, data90_09, OS_90_09, test_frequencies,
                         'tab:blue','mediumblue',
                         label='log$M_\mathrm{c} = 9.0$, $\zeta = 0.9$')
    l90_10 = plotaviolin(axs_90, data90_10, OS_90_10, test_frequencies,
                         'navy','midnightblue',
                         label='log$M_\mathrm{c} = 9.0$, $\zeta = 1.0$')
    
    plotconfig(fig_90, axs_90, fgw)
    plt.suptitle('log Mc = 9.0', fontsize=20)
    plt.legend([l90_08[0], l90_09[0], l90_10[0]],
               [l90_08[1], l90_09[1], l90_10[1]],
               bbox_to_anchor = (1.5,1.))
    plt.savefig(outdir_figure + 'logMc9.0.png', dpi=400, bbox_inches='tight')
    plt.show()
    #--------------------------------------------------------------------------
    
    #-- logMc = 8.5 -----------------------------------------------------------
    fig_85, axs_85 = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
    #plt.subplots_adjust(wspace=0.05)
      
    
    l85_08 = plotaviolin(axs_85, data85_08, OS_85_08, test_frequencies,
                         'plum','mediumorchid',
                         label='log$M_\mathrm{c} = 8.5$, $\zeta = 0.8$')
    l85_09 = plotaviolin(axs_85, data85_09, OS_85_09, test_frequencies,
                                      'mediumorchid','darkviolet',
                                      label='log$M_\mathrm{c} = 8.5$, $\zeta = 0.9$')
    l85_10 = plotaviolin(axs_85, data85_10, OS_85_10, test_frequencies,
                                      'mediumpurple','rebeccapurple',
                                      label='log$M_\mathrm{c} = 8.5$, $\zeta = 1.0$')
    
    plotconfig(fig_85, axs_85, fgw)
    
    plt.suptitle('log Mc = 8.5', fontsize=20)
    plt.legend([l85_08[0], l85_09[0], l85_10[0]],
               [l85_08[1], l85_09[1], l85_10[1]],
               bbox_to_anchor = (1.5,1.))
    plt.savefig(outdir_figure + 'logMc8.5.png', dpi=400, bbox_inches='tight')
    plt.show()
    #--------------------------------------------------------------------------
    
if __name__ == '__main__':
    main()