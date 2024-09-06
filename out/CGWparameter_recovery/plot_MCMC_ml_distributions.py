# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 16:19:29 2024

@author: kgrun
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import scipy.constants as sc

SOLAR2S = sc.G / sc.c**3 *1.98855e30
MPC2S = sc.parsec / sc.c * 1e6

def injected_values(fgw, logMc):
    Mc = (10**logMc)*SOLAR2S
    dlum = 15*MPC2S
    h = 2 * Mc**(5/3) * (np.pi*fgw)**(2/3) /dlum
    return {'cos_gwtheta': np.cos(np.pi/2),
            'gwphi': np.pi,
            'log10_h': np.log10(h),
            'log10_fgw': np.log10(fgw),
            'log10_Mc': logMc,
            'cos_inc': np.cos(np.pi),
            'phase0': np.pi,
            'psi': np.pi
            }

plt.rcParams['font.family'] = "serif"
plt.rcParams['font.sans-serif'] = "Times"

plt.rcParams['text.usetex']= False
plt.rcParams['xtick.labelsize'] = 11.0
plt.rcParams['ytick.labelsize'] = 11.0
plt.rcParams['axes.labelsize'] = 14.0




colors_85= ['plum', 'mediumorchid', 'mediumpurple', 'pink']  # 0.8, 0.9, 1.0
colors_90= ['skyblue', 'tab:blue', 'navy']  # 0.8, 0.9, 1.0
colors_95= ['lawngreen', 'limegreen', 'green']  # 0.8, 0.9, 1.0

colors_85= ['mediumpurple', 'mediumorchid', 'plum', 'pink']  # 0.8, 0.9, 1.0
colors_86= ['gray', 'lightgray']
colors_88= ['darkred', 'firebrick', 'indianred', 'lightcoral']
colors_90= ['navy', 'tab:blue', 'skyblue', 'lightskyblue']  # 0.8, 0.9, 1.0
colors_95= ['green', 'limegreen', 'lawngreen', 'lightgreen']  # 0.8, 0.9, 1.0

color_tabl = [colors_85, colors_90, colors_95]


ylabels = {'cos_gwtheta': r'cos $\theta_\mathrm{gw}$',
           'gwphi': '$\phi_\mathrm{gw}$',
           'log10_h': 'log$_{10}h$',
           'log10_fgw': 'log$_{10}f_\mathrm{gw}$',
           'log10_Mc': 'log$_{10}M_\mathrm{c}$',
           'cos_inc': r'cos $i$',
           'phase0': r'$\Phi_0$',
           'psi': '$\Psi$'
           }


basepath ='./out/CGWparameter_recovery/isotropic/normal/'


parameter = 'cos_gwtheta'


### 8.0 9.0 9.5 no pulsar term recovery ###
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
plt.subplots_adjust(wspace=0.35)
for i, lmc in enumerate([8.5, 9.0, 9.5]):
    ivs = injected_values(22.3e-9, lmc)
    
    legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_tabl[i]]
    
    for j, pd in enumerate([1.0, 1.5, 2.0]):
        data = np.genfromtxt(basepath + 'parameters_noPT_lmc{}_pd{}.txt'.format(lmc,pd), names=True)
       
        axs[i].errorbar(pd*np.ones(len(data[parameter])), data[parameter], ls='', marker='x', color='k')        #yerr=[data[parameter]-data[parameter+'_ll'], data[parameter+'_ul']-data[parameter]]
        violin = axs[i].violinplot(data[parameter].transpose(), [pd], showextrema=False, widths=0.2)
        for vl in violin['bodies']:
            vl.set_facecolor(color_tabl[i][j])
            vl.set_edgecolor(color_tabl[i][j])
            vl.set_alpha(0.5)
        
        uncertainty = data[parameter+'_ul'] - data[parameter+'_ll']
        print(lmc, pd, np.nanmean(uncertainty, axis=0), np.count_nonzero(~np.isnan(uncertainty)))

    print()
    axs[i].axhline(ivs[parameter], ls=':', color='gray')
    
    
    axs[i].set_title('log$_{10}M_\mathrm{c}$ = ' + '{}'.format(lmc))
    axs[i].set_xlabel('pulsar distance / kpc')
    axs[i].set_ylabel(ylabels[parameter])

plt.suptitle('No Pulsar Term\n', fontsize=16, y=1.003)
#plt.savefig(basepath + '/figures/noPT_{}.png'.format(parameter), bbox_inches='tight')
plt.show()



'''
lmc=8.8
ivs = injected_values(22.3e-9, lmc)
#legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_8
for j, pd in enumerate([1.0, 1.5, 2.0, 4.0]):
    data = np.genfromtxt(basepath + 'parameters_noPT_lmc{}_pd{}.txt'.format(lmc,pd), names=True)
       
        
    plt.errorbar(pd*np.ones(len(data[parameter])), data[parameter], ls='', marker='x', color='k')        #yerr=[data[parameter]-data[parameter+'_ll'], data[parameter+'_ul']-data[parameter]]
    violin = plt.violinplot(data[parameter].transpose(), [pd], showextrema=False, widths=0.2)
    for vl in violin['bodies']:
        vl.set_facecolor(colors_88[j])
        vl.set_edgecolor(colors_88[j])
        vl.set_alpha(0.5)
        
        
    plt.axhline(ivs[parameter], ls=':', color='gray')
    
plt.xlabel('pulsar distance / kpc')
plt.ylabel(ylabels[parameter])

plt.suptitle('No Pulsar Term\n', fontsize=16, y=1.003)
#plt.savefig(basepath + '/figures/noPT_{}.png'.format(parameter), bbox_inches='tight')
plt.show()
'''


lmc=8.6
ivs = injected_values(22.3e-9, lmc)
#legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_8
for j, pd in enumerate([1.0, 4.0]):
    data = np.genfromtxt(basepath + 'parameters_noPT_lmc{}_pd{}.txt'.format(lmc,pd), names=True)
       
        
    plt.errorbar(pd*np.ones(len(data[parameter])), data[parameter], ls='', marker='x', color='k')        #yerr=[data[parameter]-data[parameter+'_ll'], data[parameter+'_ul']-data[parameter]]
    violin = plt.violinplot(data[parameter].transpose(), [pd], showextrema=False, widths=0.2)
    for vl in violin['bodies']:
        vl.set_facecolor(colors_86[j])
        vl.set_edgecolor(colors_86[j])
        vl.set_alpha(0.5)
        
        
    plt.axhline(ivs[parameter], ls=':', color='gray')
    
plt.xlabel('pulsar distance / kpc')
plt.ylabel(ylabels[parameter])

plt.suptitle('No Pulsar Term\n', fontsize=16, y=1.003)
#plt.savefig(basepath + '/figures/noPT_{}.png'.format(parameter), bbox_inches='tight')
plt.show()



'''
### plots with Pulsar Term ####################################################

# pulsar phase only 

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
plt.subplots_adjust(wspace=0.35)
for i, lmc in enumerate([8.5, 9.0, 9.5]):
    ivs = injected_values(22.3e-9, lmc)
    
    legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_tabl[i]]
    
    for j, pd in enumerate([1.0, 1.5, 2.0]):
        data = np.genfromtxt(basepath + 'parameters_PT_pphase_lmc{}_pd{}.txt'.format(lmc,pd), names=True)

        axs[i].errorbar(pd*np.ones(len(data[parameter])), data[parameter], ls='', marker='x', color='k')        #yerr=[data[parameter]-data[parameter+'_ll'], data[parameter+'_ul']-data[parameter]]
        violin = axs[i].violinplot(data[parameter].transpose(), [pd], showextrema=False, widths=0.2)
        for vl in violin['bodies']:
            vl.set_facecolor(color_tabl[i][j])
            vl.set_edgecolor(color_tabl[i][j])
            vl.set_alpha(0.5)
    
    axs[i].axhline(ivs[parameter], ls=':', color='gray')
    axs[i].set_title('log$_{10}M_\mathrm{c}$ = ' + '{}'.format(lmc))
    axs[i].set_xlabel('pulsar distance / kpc')
    axs[i].set_ylabel(ylabels[parameter])
    
plt.suptitle('Pulsar Term pphase\n', fontsize=16, y=1.003)
plt.savefig(basepath + '/figures/PT_pphase_{}.png'.format(parameter), bbox_inches='tight')
plt.show()




# pulsar distance only

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
plt.subplots_adjust(wspace=0.35)
for i, lmc in enumerate([8.5, 9.0, 9.5]):
    ivs = injected_values(22.3e-9, lmc)
    
    legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_tabl[i]]
    
    for j, pd in enumerate([1.0, 1.5, 2.0]):
        data = np.genfromtxt(basepath + 'parameters_PT_pdist_lmc{}_pd{}.txt'.format(lmc,pd), names=True)

        axs[i].errorbar(pd*np.ones(len(data[parameter])), data[parameter], ls='', marker='x', color='k')        #yerr=[data[parameter]-data[parameter+'_ll'], data[parameter+'_ul']-data[parameter]]
        violin = axs[i].violinplot(data[parameter].transpose(), [pd], showextrema=False, widths=0.2)
        for vl in violin['bodies']:
            vl.set_facecolor(color_tabl[i][j])
            vl.set_edgecolor(color_tabl[i][j])
            vl.set_alpha(0.5)
    
    axs[i].axhline(ivs[parameter], ls=':', color='gray')
    axs[i].set_title('log$_{10}M_\mathrm{c}$ = ' + '{}'.format(lmc))
    axs[i].set_xlabel('pulsar distance / kpc')
    axs[i].set_ylabel(ylabels[parameter])
    
plt.suptitle('Pulsar Term pdist\n', fontsize=16, y=1.003)
plt.savefig(basepath + '/figures/PTpdist_{}.png'.format(parameter), bbox_inches='tight')
plt.show()
'''





'''
### masked plots, individually ################################################

# no PT

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
plt.subplots_adjust(wspace=0.35)


for i, lmc in enumerate([8.5, 9.0, 9.5]):
    ivs = injected_values(22.3e-9, lmc)
    
    legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_tabl[i]]
    
    for j, pd in enumerate([1.0, 1.5, 2.0]):
        data = np.genfromtxt(basepath + 'parameters_noPT_lmc{}_pd{}.txt'.format(lmc,pd), names=True)
        #mask = np.argwhere((data['log10_fgw'] < -7.6) & (data['log10_fgw'] > -7.7))
        mask = np.argwhere((data['cos_gwtheta'] < 0.4) & (data['cos_gwtheta'] > -0.4))
        
        axs[i].errorbar(pd*np.ones(len(data[parameter][mask])), data[parameter][mask], ls='', marker='x', color='k')        #yerr=[data[parameter]-data[parameter+'_ll'], data[parameter+'_ul']-data[parameter]]
        violin = axs[i].violinplot(data[parameter][mask], [pd], showextrema=False, widths=0.2)
        for vl in violin['bodies']:
            vl.set_facecolor(color_tabl[i][j])
            vl.set_edgecolor(color_tabl[i][j])
            vl.set_alpha(0.5)
    
    axs[i].axhline(ivs[parameter], ls=':', color='gray')
    
    axs[i].set_title('log$_{10}M_\mathrm{c}$ = ' + '{}'.format(lmc))
    axs[i].set_xlabel('pulsar distance / kpc')
    axs[i].set_ylabel(ylabels[parameter])

plt.suptitle('No Pulsar Term\n', fontsize=16, y=1.003)
plt.savefig(basepath + '/figures/noPT_mask_{}.png'.format(parameter), bbox_inches='tight')
plt.show()


# PT

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
plt.subplots_adjust(wspace=0.35)
for i, lmc in enumerate([8.5, 9.0, 9.5]):
    ivs = injected_values(22.3e-9, lmc)
    
    legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_tabl[i]]
    
    for j, pd in enumerate([1.0, 1.5, 2.0]):
        data = np.genfromtxt(basepath + 'parameters_PT_new_lmc{}_pd{}.txt'.format(lmc,pd), names=True)
        mask = np.argwhere((data['log10_fgw'] < -7.6) & (data['log10_fgw'] > -7.7))
        axs[i].errorbar(pd*np.ones(len(data[parameter][mask])), data[parameter][mask], ls='', marker='x', color='k')        #yerr=[data[parameter]-data[parameter+'_ll'], data[parameter+'_ul']-data[parameter]]
        violin = axs[i].violinplot(data[parameter][mask], [pd], showextrema=False, widths=0.2)
        for vl in violin['bodies']:
            vl.set_facecolor(color_tabl[i][j])
            vl.set_edgecolor(color_tabl[i][j])
            vl.set_alpha(0.5)
    
    axs[i].axhline(ivs[parameter], ls=':', color='gray')
    axs[i].set_title('log$_{10}M_\mathrm{c}$ = ' + '{}'.format(lmc))
    axs[i].set_xlabel('pulsar distance / kpc')
    axs[i].set_ylabel(ylabels[parameter])
    
plt.suptitle('Pulsar Term\n', fontsize=16, y=1.003)
plt.savefig(basepath + '/figures/PT_mask{}.png'.format(parameter), bbox_inches='tight')
plt.show()
'''








'''
### pulsar term and no pulsar term together ###################################

fgw = 22.3e-9
parameter = 'log10_fgw'

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15, 8), sharex=True)

plt.subplots_adjust(wspace=0.35, hspace=0.05)
    
for i, lmc in enumerate([8.5, 9.0, 9.5]):
    ivs = injected_values(fgw, lmc)
    
    legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_tabl[i]]
    
    for j, pd in enumerate([1.0, 1.5, 2.0]):
        data = np.genfromtxt(basepath + 'parameters_noPT_lmc{}_pd{}.txt'.format(lmc,pd), names=True)
        
        axs[0][i].errorbar(pd*np.ones(len(data[parameter])), data[parameter], ls='', marker='x', color='k')        #yerr=[data[parameter]-data[parameter+'_ll'], data[parameter+'_ul']-data[parameter]]
        violin = axs[0][i].violinplot(data[parameter].transpose(), [pd], showextrema=False, widths=0.2)
        
        #mask = np.argwhere((data['log10_fgw'] < -7.6) & (data['log10_fgw'] > -7.7))
        #axs[0][i].errorbar(pd*np.ones(len(data[parameter][mask])), data[parameter][mask], ls='', marker='x', color='k')        #yerr=[data[parameter]-data[parameter+'_ll'], data[parameter+'_ul']-data[parameter]]
        #violin = axs[0][i].violinplot(data[parameter][mask], [pd], showextrema=False, widths=0.2)
       
        for vl in violin['bodies']:
            vl.set_facecolor(color_tabl[i][j])
            vl.set_edgecolor(color_tabl[i][j])
            vl.set_alpha(0.5)
        
        
    axs[0][i].axhline(ivs[parameter], ls=':', color='gray')    
    axs[0][i].set_title('log$_{10}M_\mathrm{c}$ = ' + '{}'.format(lmc))
    
axs[0][0].set_ylabel(ylabels[parameter])
#axs[0][0].text(x=0.2, y=-14.3, s='w/o Pulsar Term', rotation='vertical', fontsize=18)

for i, lmc in enumerate([8.5, 9.0, 9.5]):
    ivs = injected_values(fgw, lmc)
    
    legend_colors = [mpatches.Patch(facecolor=cl, edgecolor=cl, alpha=0.5) for cl in color_tabl[i]]
    
    for j, pd in enumerate([1.0, 1.5, 2.0]):
        data = np.genfromtxt(basepath + 'parameters_PT_lmc{}_pd{}.txt'.format(lmc,pd), names=True)

        axs[1][i].errorbar(pd*np.ones(len(data[parameter])), data[parameter], ls='', marker='x', color='k')        #yerr=[data[parameter]-data[parameter+'_ll'], data[parameter+'_ul']-data[parameter]]
        violin = axs[1][i].violinplot(data[parameter].transpose(), [pd], showextrema=False, widths=0.2)
        
        #mask = np.argwhere((data['log10_fgw'] < -7.64) & (data['log10_fgw'] > -7.67))
        #axs[1][i].errorbar(pd*np.ones(len(data[parameter][mask])), data[parameter][mask], ls='', marker='x', color='k')        #yerr=[data[parameter]-data[parameter+'_ll'], data[parameter+'_ul']-data[parameter]]
        #violin = axs[1][i].violinplot(data[parameter][mask], [pd], showextrema=False, widths=0.2)
       
        for vl in violin['bodies']:
            vl.set_facecolor(color_tabl[i][j])
            vl.set_edgecolor(color_tabl[i][j])
            vl.set_alpha(0.5)
    
    axs[1][i].axhline(ivs[parameter], ls=':', color='gray')
    axs[1][i].set_xlabel('pulsar distance / kpc')

axs[1][0].set_ylabel(ylabels[parameter])
#axs[1][0].text(x=0.2, y=-14.25, s='w/ Pulsar Term', rotation='vertical', fontsize=18)
    
#plt.suptitle('Pulsar Term pphase\n', fontsize=16, y=1.003)
#plt.savefig(basepath + '/figures/PT_pphase_{}.png'.format(parameter), bbox_inches='tight')
plt.show()

'''