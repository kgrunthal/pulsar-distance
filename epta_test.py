# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 12:51:42 2024

@author: kgrun
"""


import numpy as np
import matplotlib.pyplot as plt
import libstempo as T
import libstempo.plot as LP, libstempo.toasim as LT
import glob
import math
import json
from enterprise.signals import gp_signals
from enterprise_extensions import model_utils, blocks
import dynesty
from enterprise.signals import signal_base
from enterprise.pulsar import Pulsar
from enterprise_extensions.frequentist import optimal_statistic as opt_stat
import enterprise.signals.parameter as parameter
from enterprise.signals import white_signals




parfolder = '/hercules/u/kgrunthal/HD/par/EPTA/'
outfolder = '/hercules/u/kgrunthal/HD/MCMCout_EPTA/'

parfiles = sorted(glob.glob(parfolder + '/*.par'))
Npsr = len(parfiles)
print(parfiles)



psrs = []
Amp = 2e-15
gamma = 13./3.

for ii in range(0,Npsr):

    # years of observations
    psr = LT.fakepulsar(parfile=parfiles[ii],
            obstimes=np.arange(53000,53000+10*365.25,28.), toaerr=0.1)

    # We now remove the computed residuals from the TOAs, obtaining (in effect) a perfect realization of the deterministic timing model. The pulsar parameters will have changed somewhat, so `make_ideal` calls `fit()` on the pulsar object.
    LT.make_ideal(psr)

    #Generate white noise
    LT.add_efac(psr,efac=1.0)

    # add to list
    psrs.append(psr)

#LT.createGWB(psrs, Amp=Amp, gam=gamma, corr='HD')
LT.createGWB(psrs, Amp=Amp, gam=gamma, turnover=True, beta = 5/3, f0 = 10**(-7.9), corr='HD')

for Psr in psrs:

    Psr.savepar(outfolder + Psr.name + '.par')
    Psr.savetim(outfolder + Psr.name + '.tim')
    T.purgetim(outfolder + Psr.name + '.tim')
    
    
    
parfiles = sorted(glob.glob(outfolder + '/*.par'))
timfiles = sorted(glob.glob(outfolder + '/*.tim'))

psrs = []
ephemeris = None
for p, t in zip(parfiles, timfiles):
    psr = Pulsar(p, t, ephem=ephemeris)
    psrs.append(psr)
    
    
# find the maximum time span to set GW frequency sampling
Tspan = model_utils.get_tspan(psrs)

# Here we build the signal model
# First we add the timing model
s = gp_signals.TimingModel()

# Then we add the white noise
# We use different white noise parameters for every backend/receiver combination
# The white noise parameters are held constant
efac = parameter.Constant(1.0)
s += white_signals.MeasurementNoise(efac=efac)

# Finally, we add the common red noise, which is modeled as a Fourier series with 30 frequency components
# The common red noise has a power-law PSD with spectral index of 4.33
s += blocks.common_red_noise_block(psd='powerlaw', prior='log-uniform', Tspan=Tspan,
                                   components=30, name='gw_crn', orf = None)

# We set up the PTA object using the signal we defined above and the pulsars
pta = signal_base.PTA([s(p) for p in psrs])

#correlation curve assuming the injected parameters
ml_params = {"gw_crn_gamma": gamma, "gw_crn_log10_A": np.log10(Amp)}
ostat = opt_stat.OptimalStatistic(psrs, pta=pta, orf='hd')
xi, rho, sig, OS, OS_sig = ostat.compute_os(params=ml_params)

print(OS, OS_sig)

# Plot the cross-correlations and compare to the Hellings-Downs curve
# Before plotting, we need to bin the cross-correlations

def weightedavg(rho, sig):
    weights, avg = 0., 0.
    for r,s in zip(rho,sig):
        weights += 1./(s*s)
        avg += r/(s*s)

    return avg/weights, np.sqrt(1./weights)

def bin_crosscorr(zeta, xi, rho, sig):

    rho_avg, sig_avg = np.zeros(len(zeta)), np.zeros(len(zeta))

    for i,z in enumerate(zeta[:-1]):
        myrhos, mysigs = [], []
        for x,r,s in zip(xi,rho,sig):
            if x >= z and x < (z+10.):
                myrhos.append(r)
                mysigs.append(s)
        rho_avg[i], sig_avg[i] = weightedavg(myrhos, mysigs)

    return rho_avg, sig_avg

# sort the cross-correlations by xi
idx = np.argsort(xi)

xi_sorted = xi[idx]
rho_sorted = rho[idx]
sig_sorted = sig[idx]

# bin the cross-correlations so that there are the same number of pairs per bin
npairs = 40

xi_mean = []
xi_err = []

rho_avg = []
sig_avg = []

i = 0
while i < len(xi_sorted):

    xi_mean.append(np.mean(xi_sorted[i:npairs+i]))
    xi_err.append(np.std(xi_sorted[i:npairs+i]))

    r, s = weightedavg(rho_sorted[i:npairs+i], sig_sorted[i:npairs+i])
    rho_avg.append(r)
    sig_avg.append(s)

    i += npairs

xi_mean = np.array(xi_mean)
xi_err = np.array(xi_err)

def get_HD_curve(zeta):

    coszeta = np.cos(zeta*np.pi/180.)
    xip = (1.-coszeta) / 2.
    HD = 3.*( 1./3. + xip * ( np.log(xip) -1./6.) )

    return HD/2


plt.plot(xi*180/np.pi, rho, '.')
zeta = np.linspace(0.01,180,100)
HD = get_HD_curve(zeta+1)
plt.plot(zeta, OS*HD, ls='--', label='Hellings-Downs', color='red', lw=1.5)
plt.errorbar(xi_mean*180/np.pi, rho_avg, xerr=xi_err*180/np.pi, color = "lightblue", lw = 5)

plt.savefig(outfolder+'HD.png')
plt.clf()


