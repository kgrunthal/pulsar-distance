# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 10:53:50 2024

@author: kgrun
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np

def plotaviolin(axs, data, frequencies, datatype='SN'):
    SN = []
    OS = []
    dOS = []
    
    OS_err = []
    OS_max, OS_range = [], []
    
    for i, f in enumerate(frequencies):
        OS.append(data[int(3*i)])
        SN.append(data[int(3*i+2)])
        #SN.append(data[int(3*i)]/np.std(data[int(3*i)]))
        
        
        hist_OS, binedge_OS = np.histogram(data[int(3*i)], bins='auto')
        loc_OS = (binedge_OS[np.argmax(hist_OS)]+binedge_OS[np.argmax(hist_OS)+1])/2
        OS_max.append(loc_OS)
        OS_range.append(np.std(data[int(3*i)]))
        OS_err.append(np.mean(data[int(3*i+1)]))
        
    OS = np.array(OS)
    SN = np.array(SN)
    
    
    #axs[0].errorbar(np.log10(frequencies), OS_max, yerr=OS_range, ls='', fmt='rx', capsize=2)
    #axs[0].errorbar(np.log10(frequencies)+0.02  , OS_max, yerr=OS_err, ls='', fmt='bx', capsize=2)
    #print(OS_range[6],OS_err[6])
    #print(OS_range[0],OS_err[0])
    
    
    # plot: OS values
    if datatype == 'OS':
        violin_OS = axs[1].violinplot(OS.transpose(), np.log10(frequencies),
                            showextrema=False, widths=0.05)
        for vl in violin_OS['bodies']:
            vl.set_facecolor('darkcyan')
            vl.set_edgecolor('darkcyan')
            vl.set_alpha(0.5)
    
    # right plot: SN values
    elif datatype == 'SN': 
        violin_SN = axs[1].violinplot(SN.transpose(), np.log10(frequencies),
                                      showextrema=False, widths=0.05)
        for vl in violin_SN['bodies']:
            vl.set_facecolor('darkcyan')
            vl.set_edgecolor('darkcyan')
            vl.set_alpha(0.5)
    
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
    
    free_spectrum = np.loadtxt('MCMCout_isotropic_singleRN_CGW9.5_pd1.0_1/chain_1.txt')
    OS_data_nm = np.loadtxt('out/WN_CGW_singleRN/run_1/OS_spectrum_CGW9.5_pd1.0_NM.txt')
    OS_data = np.loadtxt('out/WN_CGW_singleRN/run_1/OS_spectrum_CGW9.5_pd1.0.txt').transpose()
    
    ###########################################################################
    ###########################################################################
    
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(6, 8), sharex=True)
    plt.subplots_adjust(hspace=0.05)
    
    burn = int(0.2*len(free_spectrum[:,0]))
    chain_cut = free_spectrum[burn:,:]
    violin_fs = axs[0].violinplot(chain_cut[::32,0:20], OS_data[0],
                                  showextrema=False, widths=2e-9)
    for vl in violin_fs['bodies']:
            vl.set_facecolor('darkcyan')
            vl.set_edgecolor('darkcyan')
            vl.set_alpha(0.6)
            vl.set_zorder(2)
            
    axs[1].errorbar(OS_data[0], OS_data[1]*1e12, yerr=OS_data[2]*1e12,
                    ls='',color='darkcyan', marker='o', capsize=2, zorder=1)
    #l95_10 = plotaviolin(axs, OS_data, test_frequencies, datatype='OS')
    axs[1].set_xscale('log')
    
    
    # lines
    axs[0].axvline(22.3e-9, ls=':', color='gray', zorder=0)
    axs[1].axvline(22.3e-9, ls=':', color='gray', zorder=0)
    axs[0].text(20e-9, -7, 'CGW', rotation='vertical', color='gray')
    axs[1].text(20e-9, 1, 'CGW', rotation='vertical', color='gray')
    
    axs[0].axvline(12.77e-9, ls=':', color='gray', zorder=0)
    axs[1].axvline(12.77e-9, ls=':', color='gray', zorder=0)
    axs[0].text(11.4e-9, -7.1, 'no corr', rotation='vertical', color='gray')
    axs[1].text(11.4e-9, 0.85, 'no corr', rotation='vertical', color='gray')  
    
    
    # labels
    axs[1].set_xlabel('frequency / Hz')
    axs[0].set_ylabel(r'log$_{10}$($\rho$ / s)')
    axs[1].set_ylabel(r'OS / $\times 10^{-12}$')
    plt.savefig('comparison_freespectrum_os.pdf', bbox_inches='tight')
    plt.show()
    #--------------------------------------------------------------------------
    
    
    
    
        
if __name__ == '__main__':
    main()