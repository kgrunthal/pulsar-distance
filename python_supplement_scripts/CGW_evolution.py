# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 12:58:19 2024

@author: kgrun
"""


import argparse
import numpy as np
import pickle
import json
import glob
import os
import scipy.constants as sc
#import scipy.signal as signal
import matplotlib.pyplot as plt

from astropy import units as u

from astropy.coordinates import SkyCoord

SOLAR2S = sc.G / sc.c**3 * 1.98855e30
KPC2S = sc.parsec / sc.c * 1e3
MPC2S = sc.parsec / sc.c * 1e6


def waveform(
    TOAs,
    RA, DEC,
    gwtheta,
    gwphi,
    mc,
    dist,
    fgw,
    phase0,
    psi,
    inc,
    pdist=1.0,
    zeta = None,
    pphase=None,
    psrTerm=True,
    evolve=True,
    phase_approx=False,
    tref=0,
    ):
    """
    Function to create GW-induced residuals from a SMBMB as
    defined in Ellis et. al 2012,2013. Tries to be smart about it...

    :param psr: pulsar object
    :param gwtheta: Polar angle of GW source in celestial coords [radians]
    :param gwphi: Azimuthal angle of GW source in celestial coords [radians]
    :param mc: Chirp mass of SMBMB [solar masses]
    :param dist: Luminosity distance to SMBMB [Mpc]
    :param fgw: Frequency of GW (twice the orbital frequency) [Hz]
    :param phase0: Initial Phase of GW source [radians]
    :param psi: Polarization of GW source [radians]
    :param inc: Inclination of GW source [radians]
    :param pdist: Pulsar distance to use other than those in psr [kpc]
    :param pphase: Use pulsar phase to determine distance [radian]
    :param psrTerm: Option to include pulsar term [boolean]
    :param evolve: Option to exclude evolution [boolean]
    :param tref: Fidicuial time at which initial parameters are referenced

    :returns: Vector of induced residuals
    """

    # convert units
    mc *= SOLAR2S  # convert from solar masses to seconds
    dist *= MPC2S  # convert from Mpc to seconds

    # define initial orbital frequency
    w0 = np.pi * fgw
    phase0 /= 2  # orbital phase
    w053 = w0 ** (-5 / 3)

    # define variable for later use
    cosgwtheta, cosgwphi = np.cos(gwtheta), np.cos(gwphi)
    singwtheta, singwphi = np.sin(gwtheta), np.sin(gwphi)
    sin2psi, cos2psi = np.sin(2 * psi), np.cos(2 * psi)
    incfac1, incfac2 = 0.5 * (3 + np.cos(2 * inc)), 2 * np.cos(inc)

    # unit vectors to GW source
    m = np.array([singwphi, -cosgwphi, 0.0])
    n = np.array([-cosgwtheta * cosgwphi, -cosgwtheta * singwphi, singwtheta])
    omhat = np.array([-singwtheta * cosgwphi, -singwtheta * singwphi, -cosgwtheta])

    # various factors invloving GW parameters
    fac1 = 256 / 5 * mc ** (5 / 3) * w0 ** (8 / 3)
    fac2 = 1 / 32 / mc ** (5 / 3)
    fac3 = mc ** (5 / 3) / dist

    # pulsar location
    
    
    psr_c = SkyCoord(RA, DEC, frame='icrs')
    ptheta = np.pi/2. - psr_c.dec.value
    pphi = psr_c.ra.value
    

    # use definition from Sesana et al 2010 and Ellis et al 2012
    phat = np.array([np.sin(ptheta) * np.cos(pphi), np.sin(ptheta) * np.sin(pphi), np.cos(ptheta)])

    fplus = 0.5 * (np.dot(m, phat) ** 2 - np.dot(n, phat) ** 2) / (1 + np.dot(omhat, phat))
    fcross = (np.dot(m, phat) * np.dot(n, phat)) / (1 + np.dot(omhat, phat))
    cosMu = -np.dot(omhat, phat)

    # get values from pulsar object
    toas = TOAs * 86400 - tref
    if pphase is not None:
        pd = pphase / (2 * np.pi * fgw * (1 - cosMu)) / KPC2S
    else:
        pd = pdist

    # convert units
    pd *= KPC2S  # convert from kpc to seconds

    # get pulsar time
    tp = toas - pd * (1 - cosMu)

    # evolution
    if evolve:

        # calculate time dependent frequency at earth and pulsar
        omega = w0 * (1 - fac1 * toas) ** (-3 / 8)
        omega_p = w0 * (1 - fac1 * tp) ** (-3 / 8)        

        # calculate time dependent phase
        phase = phase0 + fac2 * (w053 - omega ** (-5 / 3))
        phase_p = phase0 + fac2 * (w053 - omega_p ** (-5 / 3))

    # use approximation that frequency does not evlolve over observation time
    elif phase_approx:

        # frequencies and phases
        omega = w0
        phase = phase0 + omega * toas
        if zeta < 0.:
            omega_p = w0 * (1 + fac1 * pd * (1 - cosMu)) ** (-3 / 8)
            phase_p = phase0 + fac2 * (w053 - omega_p ** (-5 / 3)) + omega_p * toas #+ np.random.uniform(0,2*np.pi)
            
        elif zeta == 1.:
            omega_p = omega
            phase_p = np.random.uniform(0,2*np.pi) + omega_p * toas
        else:
            omega_p = zeta*float(omega)
            phase_p = np.random.uniform(0,2*np.pi) + omega_p * toas


    # no evolution
    else:

        # monochromatic
        omega = w0
        omega_p = omega

        # phases
        phase = phase0 + omega * toas
        phase_p = phase0 + omega * tp

    # define time dependent coefficients
    At = np.sin(2 * phase) * incfac1
    Bt = np.cos(2 * phase) * incfac2
    At_p = np.sin(2 * phase_p) * incfac1
    Bt_p = np.cos(2 * phase_p) * incfac2

    # now define time dependent amplitudes
    alpha = fac3 / omega ** (1 / 3)
    alpha_p = fac3 / omega_p ** (1 / 3)

    # define rplus and rcross
    rplus = alpha * (At * cos2psi + Bt * sin2psi)
    rcross = alpha * (-At * sin2psi + Bt * cos2psi)
    rplus_p = alpha_p * (At_p * cos2psi + Bt_p * sin2psi)
    rcross_p = alpha_p * (-At_p * sin2psi + Bt_p * cos2psi)

    # residuals
    if psrTerm:
        res = fplus * (rplus_p - rplus) + fcross * (rcross_p - rcross)
        res_p = fplus * rplus_p + fcross * rcross_p
        res_e = -fplus * rplus - fcross * rcross
        
    else:
        res = -fplus * rplus - fcross * rcross


    return [res, res_e, res_p], [omega, omega_p]





# Pulsar #
name='J1514-1428'
#RA, DEC = '00h24m47s', '-20d29m14s'       #J0024-2029
RA, DEC = '15h14m47s', '-14d28m14s'       #J1514-1428


TOAs = np.arange(50000, 53652, 14)
Tobs = (np.max(TOAs) - np.min(TOAs)) * 24*3600
fobs = 1/Tobs
pdist = 1.

print(fobs)
# CGW #
gwtheta = np.pi/2.
gwphi = np.pi

mc = 10**9.5
dl = 15         #Mpc
fgw = 22.3*1e-9

phase0 = np.pi
psi = np.pi
inc = np.pi

psrTerm = True
evolve = False
phase_approx = True
zeta = -1


residuals, omega = waveform(TOAs, RA, DEC, gwtheta, gwphi, mc, dl, fgw, phase0, psi, inc,
                     pdist=1.0, zeta=zeta, psrTerm = psrTerm, evolve=evolve, phase_approx=phase_approx)

plt.plot(TOAs, residuals[0], color='k', ls = '-', lw=2)
plt.plot(TOAs, residuals[1], color='b', ls='--', label='Earth')
plt.plot(TOAs, residuals[2], color = 'purple', ls = ':', label='Pulsar')

plt.legend(loc='best')
plt.show()





### evolve = False #####


colors = ['mediumorchid', 'tab:blue', 'limegreen']  # 8.5, 9.0, 9.5
distances = np.arange(0.5, 5, 0.01)


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7, 9 y), sharex=True)
plt.subplots_adjust(hspace=0.05)

fgw_1 = 5e-9  #Hz
for ii, lmc in enumerate([8.5,9.0,9.5]):
    f_e, f_p = np.zeros(len(distances)), np.zeros(len(distances))
    
    for i, dist in enumerate(distances):        
        residuals, omega = waveform(TOAs, RA, DEC, gwtheta, gwphi, 10**lmc, dl, fgw_1, phase0, psi, inc,
                                pdist=dist, zeta=zeta, psrTerm = psrTerm, evolve=evolve, phase_approx=phase_approx)
        f_e[i] = omega[0]/np.pi
        f_p[i] = omega[1]/np.pi
        
    axs[1].plot(distances, 1e9*f_p, label = 'pulsar term, log$_{10}M_\mathrm{c}$='+' {}'.format(lmc), color=colors[ii])
axs[1].plot(distances, 1e9*f_e, label = 'earth term', color='k')
    
for f in range(1,3):
    axs[1].axhline(1e9*f*fobs, ls=':', color='gray')


fgw_2 = 22.3e-9  #Hz
for ii, lmc in enumerate([8.5,9.0,9.5]):
    f_e, f_p = np.zeros(len(distances)), np.zeros(len(distances))
    
    for i, dist in enumerate(distances):        
        residuals, omega = waveform(TOAs, RA, DEC, gwtheta, gwphi, 10**lmc, dl, fgw_2, phase0, psi, inc,
                                pdist=dist, zeta=zeta, psrTerm = psrTerm, evolve=evolve, phase_approx=phase_approx)
        f_e[i] = omega[0]/np.pi
        f_p[i] = omega[1]/np.pi
        
    axs[0].plot(distances, 1e9*f_p, label = 'pulsar term, log$_{10}M_\mathrm{c}$='+' {}'.format(lmc), color=colors[ii])
axs[0].plot(distances, 1e9*f_e, label = 'earth term', color='k')
    
for f in range(2,8):
    axs[0].axhline(1e9*f*fobs, ls=':', color='gray')


axs[0].set_xlim(1.5, 5)
axs[1].set_xlim(1.5, 5)
axs[1].set_xlabel('Pulsar distance / kpc')
axs[0].set_ylabel('$f_\mathrm{gw}$ / nHz')
axs[1].set_ylabel('$f_\mathrm{gw}$ / nHz')
plt.legend(bbox_to_anchor=(0.02,0.6))
plt.savefig('CGW_FreqEvo_{}.png'.format(name), bbox_inches ='tight', dpi=400)
plt.show()
'''


### evolve = True ####

residuals, omega = waveform(TOAs, RA, DEC, gwtheta, gwphi, mc, dl, fgw, phase0, psi, inc,
                     pdist=1.0, zeta=zeta, psrTerm = psrTerm, evolve=True, phase_approx=False)
print(min(omega[0]/np.pi), max(omega[0]/np.pi))
plt.plot(TOAs, omega[0]/np.pi, color='b', label='Earth')
plt.plot(TOAs, omega[1]/np.pi, color = 'purple', label='Pulsar')
plt.legend(loc='best')
plt.show()

'''