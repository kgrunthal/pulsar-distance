# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 11:05:26 2024

@author: kgrun
"""


import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches




def S(frequencies, log10A, gamma):
    A = 10**log10A
    fyr = 1/(365.25*24*3600)
    return A**2/(12*np.pi**2) * (frequencies/fyr)**(-1*gamma) * fyr**(-3)




def plotaviolin(data, OS_data, frequencies, facecolor, edgecolor, label=''):
    SN = []
    OS = []
    for i, f in enumerate(frequencies):
        SN.append(data[int(3*i+2)])
        OS.append(data[int(3*i)])
    SN = np.array(SN)  
    OS = np.array(OS)
    violin = plt.violinplot(SN.transpose(), np.log10(frequencies),
                            showextrema=False, widths=0.05)
    
    #plt.errorbar(np.log10(frequencies), OS_data[1]/OS_data[2],
    #             marker='o', markersize=4, ls='', color=edgecolor)
    for vl in violin['bodies']:
        vl.set_facecolor(facecolor)
        vl.set_edgecolor(edgecolor)
        vl.set_alpha(0.5)
    legend = [mpatches.Patch(facecolor=facecolor, edgecolor=edgecolor, alpha=0.5), label]
    
    return violin, legend



def plot_comparison(freespectrum, plparams, OS_fs_NM, OS_pl_NM, OS_fsfull_NM, frequencies):
    # get free spectrum values
    fs_val = np.fromiter(freespectrum.values(), dtype=float)
    pl_val = np.fromiter(plparams.values(), dtype=float)
    
    # get OS values from free spectrum analysis
    fs_SN, fs_OS  = [], []
    for i, f in enumerate(frequencies):
        fs_SN.append(OS_fs_NM[int(3*i+2)])
        fs_OS.append(OS_fs_NM[int(3*i)])
    fs_SN = np.array(fs_SN)  
    fs_OS = np.array(fs_OS)
    
    

    fig, axs = plt.subplots(nrows = 1, ncols = 4, figsize=(22, 4))
    plt.subplots_adjust(wspace=0.4)
    
    # free spectrum rhos
    axs[0].errorbar(np.log10(frequencies), fs_val, ls ='', fmt='kx')
    #axs[0].plot(np.log10(frequencies), np.log10(S(frequencies, pl_val[1], pl_val[0])))
    
    # frequency resolved free spectrum OS
    violin_fs_OS = axs[1].violinplot(np.log10(fs_OS.transpose()), np.log10(frequencies),
                                 showextrema=False, widths=0.05)
    for vl in violin_fs_OS['bodies']:
        vl.set_facecolor('lightblue')
        vl.set_edgecolor('tab:blue')
        vl.set_alpha(0.5)
    
    violin_fs_SN = axs[2].violinplot(fs_SN.transpose(), np.log10(frequencies),
                                 showextrema=False, widths=0.05)
    for vl in violin_fs_SN['bodies']:
        vl.set_facecolor('lightblue')
        vl.set_edgecolor('tab:blue')
        vl.set_alpha(0.5)

    violin_OS = axs[3].violinplot(np.array([OS_pl_NM[2], OS_fsfull_NM[2]]).transpose() , [0., 1.],
                                     showextrema=False, widths=0.05)
    for vl in violin_OS['bodies']:
        vl.set_facecolor('green')
        vl.set_edgecolor('darkgreen')
        vl.set_alpha(0.5)
    
    plot_config(fig, axs)
    return None





def distribution_spectrum(data, frequencies):
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


def distribution_single(data):
    hist, binedge = np.histogram(data, bins='auto')
    loc = (binedge[np.argmax(hist)]+binedge[np.argmax(hist)+1])/2

    return loc



def plot_config(fig, axs):
    
    axs[0].set_xlabel('log$_{10}$($f$ / Hz)' , fontsize=14)
    axs[1].set_xlabel('log$_{10}$($f$ / Hz)' , fontsize=14)
    axs[2].set_xlabel('log$_{10}$($f$ / Hz)' , fontsize=14)
    
    axs[0].set_ylabel(r'log$_{10}\rho$', fontsize=14)
    axs[1].set_ylabel('$\hat{A}_\mathrm{GW}^2$', fontsize=14)
    axs[2].set_ylabel('S/N', fontsize=14)
    axs[3].set_ylabel('S/N', fontsize=14)
    
    axs[3].set_xticks([0., 1])
    axs[3].set_xticklabels(['powerlaw', 'free spectrum'], rotation=45)
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
    
    basefolder = './out/powerlaw_vs_spectrum/last_bin/'
    GWB_type = 'GWBsinglebin_lastbin'
    realisations = 5
    
    maxlike_fs = np.zeros((realisations, bins))
    OS_spec_distr = np.zeros((realisations, bins))
    SN_spec_distr = np.zeros((realisations, bins))
    SN_fullspec_distr = np.zeros(realisations)
    SN_powerlaw_distr = np.zeros(realisations)
    
    for i in range(5):
        maxlike_out = './MCMCout_{}_{}/freespectrum/maxlike.json'.format(GWB_type, i+1)
        maxlike_freespectrum = json.load(open(maxlike_out))
        maxlike_fs[i] = np.fromiter(maxlike_freespectrum.values(), dtype=float)
        
        OS_spectrum = np.loadtxt(basefolder + 'OS_spectrum_{}_{}_NM.txt'.format(GWB_type, i+1))
        OS_fullspectrum = np.loadtxt(basefolder + 'OS_fullspectrum_{}_{}_NM.txt'.format(GWB_type, i+1)).transpose()
        OS_powerlaw = np.loadtxt(basefolder + 'OS_powerlaw_{}_{}_NM.txt'.format(GWB_type, i+1)).transpose()
        
        OS_spec_distr[i], SN_spec_distr[i] = distribution_spectrum(OS_spectrum, test_frequencies)
        
        SN_fullspec_distr[i] = distribution_single(OS_fullspectrum[2])
        SN_powerlaw_distr[i] = distribution_single(OS_powerlaw[2])
        
        
        
        
    # produce ranges for plotting
    maxlike_avg = np.average(maxlike_fs, axis=0)
    maxlike_range = np.max(maxlike_fs, axis=0) - maxlike_avg
    
    OS_fs_avg = np.average(OS_spec_distr, axis=0)
    OS_fs_range = np.max(OS_spec_distr, axis=0) - OS_fs_avg
    SN_fs_avg = np.average(SN_spec_distr, axis=0)
    SN_fs_range = np.max(SN_spec_distr, axis=0) - SN_fs_avg
    
    SN_fsfull_avg = np.average(SN_fullspec_distr, axis=0)
    SN_fsfull_range = np.max(SN_fullspec_distr, axis=0) - SN_fsfull_avg
    SN_pl_avg = np.average(SN_powerlaw_distr, axis=0)
    SN_pl_range = np.max(SN_powerlaw_distr, axis=0) - SN_pl_avg
    
    
    
    fig, axs = plt.subplots(nrows = 1, ncols = 4, figsize=(22, 4))
    plt.subplots_adjust(wspace=0.4)
    
    # free spectrum rhos
    axs[0].errorbar(np.log10(test_frequencies), maxlike_avg, yerr=maxlike_range,
                    markersize=0, capsize=3, linewidth=3, elinewidth=3, ls ='')
    #axs[0].plot(np.log10(frequencies), np.log10(S(frequencies, pl_val[1], pl_val[0])))
    
    # frequency resolved free spectrum OS
    axs[1].errorbar(np.log10(test_frequencies), np.log10(OS_fs_avg), yerr=OS_fs_range,
                    markersize=0, capsize=3, linewidth=3, elinewidth=3, ls ='')
    
    # frequency resolved free spectrum SN
    axs[2].errorbar(np.log10(test_frequencies), SN_fs_avg, yerr=SN_fs_range,
                    markersize=0, capsize=3, linewidth=3, elinewidth=3, ls ='')
    
    axs[3].errorbar([0, 1], [SN_pl_avg, SN_fsfull_avg], yerr=[SN_pl_range, SN_fsfull_range],
                    markersize=0, capsize=3, linewidth=3, elinewidth=3, ls ='')

    
    plot_config(fig, axs)
    plt.suptitle('GWB single bin')

    plt.savefig(basefolder + 'realisations_{}.png'.format(GWB_type), dpi=400, bbox_inches= 'tight')
    plt.show()
    return None
    
        
        

    
    
    
    
    
    
if __name__ == '__main__':
    main()