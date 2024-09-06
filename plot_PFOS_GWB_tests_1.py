# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 16:03:07 2024

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
    
    basefolder = './out/powerlaw_vs_spectrum/'
    run_pl = 5      # powerlaw
    run_bpl = 1     # broken powerlaw
    run_sgl = 1     # single bin
    
    ###########################################################################
    ###  noisemarginalised single OS value  ###################################
    
    GWBpl_powerlaw_NM = np.loadtxt(basefolder + 'OS_powerlaw_GWB_{}_NM.txt'.format(run_pl)).transpose()
    GWBbroken_powerlaw_NM = np.loadtxt(basefolder + 'OS_powerlaw_GWBbroken_{}_NM.txt'.format(run_bpl)).transpose()
    GWBsinglebin_powerlaw_NM = np.loadtxt(basefolder + 'OS_powerlaw_GWBsinglebin_{}_NM.txt'.format(run_sgl)).transpose()
        
    GWBpl_fullspectrum_NM = np.loadtxt(basefolder + 'OS_fullspectrum_GWB_{}_NM.txt'.format(run_pl)).transpose()
    GWBbroken_fullspectrum_NM = np.loadtxt(basefolder + 'OS_fullspectrum_GWBbroken_{}_NM.txt'.format(run_bpl)).transpose()
    GWBsinglebin_fullspectrum_NM = np.loadtxt(basefolder + 'OS_fullspectrum_GWBsinglebin_{}_NM.txt'.format(run_sgl)).transpose()
    
    
    ###########################################################################
    ###  noisemarginalised frequency OS values  ###############################
    
    GWBpl_spectrum_NM = np.loadtxt(basefolder + 'OS_spectrum_GWB_{}_NM.txt'.format(run_pl))
    GWBbroken_spectrum_NM = np.loadtxt(basefolder + 'OS_spectrum_GWBbroken_{}_NM.txt'.format(run_bpl))
    GWBsinglebin_spectrum_NM = np.loadtxt(basefolder + 'OS_spectrum_GWBsinglebin_{}_NM.txt'.format(run_sgl))
    
    ###########################################################################
    ###  maxlike values  ######################################################
    maxlike_freespectrum_GWBpl = json.load(open('./MCMCout_GWB_{}/freespectrum/maxlike.json'.format(run_pl)))
    maxlike_freespectrum_GWBbroken = json.load(open('./MCMCout_GWBbroken_{}/freespectrum/maxlike.json'.format(run_bpl)))
    maxlike_freespectrum_GWBsinglebin = json.load(open('./MCMCout_GWBsinglebin_{}/freespectrum/maxlike.json'.format(run_sgl)))
    
    maxlike_powerlaw_GWBpl = json.load(open('./MCMCout_GWB_{}/powerlaw/maxlike.json'.format(run_pl)))
    maxlike_powerlaw_GWBbroken = json.load(open('./MCMCout_GWBbroken_{}/powerlaw/maxlike.json'.format(run_bpl)))
    maxlike_powerlaw_GWBsinglebin = json.load(open('./MCMCout_GWBsinglebin_{}/powerlaw/maxlike.json'.format(run_sgl)))
    
    ###########################################################################
    ###########################################################################
    outdir = './out/powerlaw_vs_spectrum/'
    plot_comparison(maxlike_freespectrum_GWBpl, maxlike_powerlaw_GWBpl,
                    GWBpl_spectrum_NM,
                    GWBpl_powerlaw_NM, GWBpl_fullspectrum_NM,
                    test_frequencies)
    plt.suptitle('Powerlaw GWB')
    plt.savefig(outdir + 'gwb_powerlaw_{}.png'.format(run_pl), bbox_inches='tight', dpi=400)
    plt.show()
    
    #--------------------------------------------------------------------------
    plot_comparison(maxlike_freespectrum_GWBbroken, maxlike_powerlaw_GWBbroken,
                    GWBbroken_spectrum_NM,
                    GWBbroken_powerlaw_NM, GWBbroken_fullspectrum_NM,
                    test_frequencies)
    plt.suptitle('Broken powerlaw GWB')
    plt.savefig(outdir + 'gwb_broken_powerlaw_{}.png'.format(run_bpl), bbox_inches='tight', dpi=400)
    plt.show()
    
    #--------------------------------------------------------------------------
    plot_comparison(maxlike_freespectrum_GWBsinglebin, maxlike_powerlaw_GWBsinglebin,
                    GWBsinglebin_spectrum_NM,
                    GWBsinglebin_powerlaw_NM, GWBsinglebin_fullspectrum_NM,
                    test_frequencies)
    plt.suptitle('Single bin GWB')
    plt.savefig(outdir + 'gwb_singlebin_{}.png'.format(run_sgl), bbox_inches='tight', dpi=400)
    plt.show()
    
    
    
    
    
if __name__ == '__main__':
    main()