#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 12:09:11 2023

@author: kgrunthal
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

import ephem

import libstempo as T
import libstempo.toasim as LT
import libstempo.plot as LP

from enterprise.signals import signal_base, gp_signals, white_signals, utils, gp_priors
from enterprise.pulsar import Pulsar
import enterprise.signals.parameter as parameter

from enterprise_extensions import model_utils, blocks, deterministic
from enterprise_extensions.frequentist import optimal_statistic as opt_stat
from enterprise_extensions import sampler as sp

import corner
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc

import h5py


SOLAR2S = sc.G / sc.c**3 * 1.98855e30
KPC2S = sc.parsec / sc.c * 1e3
MPC2S = sc.parsec / sc.c * 1e6



def S(A, gamma, f):
    fc = 1/(365.25*86000)
    A_c = (A**2)/(12*np.pi**2.)
    return A_c*((f/fc)**(-1*gamma))*(fc**(-3.))

def singlebin(amp, gamma, f_window, f):
    spectrum = 1e-50*np.ones(len(f))
    window = [ii for ii in range(len(f)) if f[ii] > f_window[0] and f[ii] < f_window[1]]
    spectrum[window] = S(amp, gamma, f[window])
    return spectrum




def ADD_CGW(
    psr,
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
    if "RAJ" and "DECJ" in psr.pars(which='all'):
        ptheta = np.pi / 2 - psr["DECJ"].val
        pphi = psr["RAJ"].val
    elif "ELONG" and "ELAT" in psr.pars():
        fac = 180.0 / np.pi
        coords = ephem.Equatorial(ephem.Ecliptic(str(psr["ELONG"].val * fac), str(psr["ELAT"].val * fac)))

        ptheta = np.pi / 2 - float(repr(coords.dec))
        pphi = float(repr(coords.ra))

    # use definition from Sesana et al 2010 and Ellis et al 2012
    phat = np.array([np.sin(ptheta) * np.cos(pphi), np.sin(ptheta) * np.sin(pphi), np.cos(ptheta)])

    fplus = 0.5 * (np.dot(m, phat) ** 2 - np.dot(n, phat) ** 2) / (1 + np.dot(omhat, phat))
    fcross = (np.dot(m, phat) * np.dot(n, phat)) / (1 + np.dot(omhat, phat))
    cosMu = -np.dot(omhat, phat)

    # get values from pulsar object
    toas = psr.toas() * 86400 - tref
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

        # frequencies
        omega = w0
        
        if zeta == None:
            omega_p = w0 * (1 + fac1 * pd * (1 - cosMu)) ** (-3 / 8)
            
        elif zeta == 1.:
            omega_p = omega
        else:
            omega_p = zeta*float(omega)

        # phases
        phase = phase0 + omega * toas
        #phase_p = phase0 + fac2 * (w053 - omega_p ** (-5 / 3)) + omega_p * toas
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
        #res = fplus * rplus_p + fcross * rcross_p  # pulsar term only
        
    else:
        res = -fplus * rplus - fcross * rcross

    psr.stoas[:] += res / 86400





def prime_pulsars(psrs, pdistances, signal, Ncgw, fgw, Mc, zeta=None, psrTerm=True, phase_approx=True, evolve=False, Filter=False, f0=None):
    ePSRs = []
    sgn = signal.split(',')
    print('psrTerm=', psrTerm)
    for ii, ltp in enumerate(psrs):
        
        LT.make_ideal(ltp)
        stoas_temp = ltp.stoas[:]
        
        outstr = ltp.name
        LT.add_efac(ltp, efac=1.0)
        outstr += ' added WN'
        
        for s in sgn:
            if s=='RN':
                log10A = np.random.uniform(-14.,-12.)
                gamma = np.random.uniform(1.7,2.2)
                LT.add_rednoise(ltp,10**log10A, gamma, components=20)
                outstr += ' added RN (logA = {}, gamma={})'.format(log10A, gamma)
            
            elif s=='CRN':
                tspan = 5.*365.25    # days
                LT.add_rednoise(ltp,10**-14,2.2, tspan=tspan, components=3)
                outstr += ' added CRN'
                
            elif s=='CGW':
                print(fgw)
                for n in range(int(Ncgw)):
                    ADD_CGW(ltp,
                            np.pi, 0., #CGWtheta CGWphi
                            Mc[n], 15., fgw[n],
                            0., 0., 0.,
                            pdist = pdistances[ii],
                            zeta = zeta,
                            psrTerm=psrTerm, phase_approx=phase_approx, evolve=evolve)
                    outstr += ' added CGW {}'.format(n+1)
            
            else:
                print('unsupported signal type')
                break
            
            ''' old methods
            elif s=='GWB':
                # old method
                gamma = 13./3.
                alpha = (3.-gamma)/2.
                log10_A_gw = -14.5
                
                # single bin
                #flow = 7e-9
                #fhigh = 1.1e-8

                # first bin
                flow = 5e-10
                fhigh = 1e-5
                LT.add_gwb(ltp, gwAmp=10**log10_A_gw, alpha=alpha, flow=flow, fhigh=fhigh)
                outstr += ' added GWB with log10A={}, alpha={:.2}'.format(log10_A_gw,alpha)  
            
            elif s=='GWBbroken':
                # old method
                gamma_1 = 16./3.
                alpha_1 = (3.-gamma_1)/2.
                log10_A_gw_1 = -14.3
                flow_1 = 1e-9
                fhigh_1 = 9.6*1e-9
                
                LT.add_gwb(ltp, gwAmp=10**log10_A_gw_1, alpha=alpha_1, flow=flow_1, fhigh=fhigh_1)
                
                gamma_2 = 10./3.
                alpha_2 = (3.-gamma_1)/2.
                log10_A_gw_2 = -14.3
                flow_2 = 9.6*1e-9
                fhigh_2 = 1e-5

                LT.add_gwb(ltp, gwAmp=10**log10_A_gw_2, alpha=alpha_2, flow=flow_2, fhigh=fhigh_2)

                outstr += ' added broken GWB '
            '''
        print(outstr)
    
    cstring = ''
    for s in sgn:
        if s=='GWB':
            amp = 2e-15
            gamma = 13./3.
                
            LT.createGWB(psrs, amp, gamma)
            cstring += ' added GWB with log10A={}, gamma={:.2}'.format(np.log10(amp),gamma)


        elif s=='GWBbroken':
            amp = 2e-15
            gamma = 13./3.
            beta = 5/3
            f0 = 10**(-7.9)
                
            LT.createGWB(psrs, amp, gamma, turnover=True,  beta = beta, f0 = f0)
            cstring += ' added broken GWB '


        elif s=='GWBsinglebin':
            amp = 2e-15
            gamma = 13./3.
            window = [3./(3625*86400), 4./(3625*86400)]
            
            f_test = np.linspace(1/(3625*86400), 20/(3625*86400), 100) 
            S = singlebin(amp, gamma, window, f_test)
            f_test = f_test.astype('float128')
            S = S.astype('float128')
                
            customSpec = np.array([f_test, S]).transpose()
            LT.createGWB(psrs, amp, gamma, userSpec = customSpec)
            cstring += ' added single bin GWB '
            
    print(cstring, flush=True)   
        
    for p in psrs:
        ePSRs.append(Pulsar(p, dist=pdistances[ii]))
        
    return ePSRs



def set_up_global_options():
    parser = argparse.ArgumentParser(description='Generate sensitivity curves for PTA datasets.')
    parser.add_argument('--par', type=str, default=None, help='Path to the parfiles.')
    parser.add_argument('--outdir', type=str, default=None, help='Directory for the MCMC output')
    parser.add_argument('--result', type=str, default=None, help='Directory and filename for the OS stats output')

    parser.add_argument('--signal', type=str, default='CGW', help='Which signals to inject to the pulsar objects')
    parser.add_argument('--lmc', type=float, nargs='*', default=9.5)
    parser.add_argument('--fgw', type=float, nargs='*', default=22.3)
    parser.add_argument('--ncgw', type=int, default=1)

    parser.add_argument('--psrTerm', action="store_true", default=False, help='Simulate CGW with or without pulsar term')
    parser.add_argument('--zeta', type=float, default=1.0)
    return parser.parse_args()






    
    
def main():
    args = set_up_global_options()
    
    ### load pulsars ###

    parfiles = glob.glob(args.par+'/*.par')

    obstimes = np.arange(50000,53652,14)
    toaerr = 0.1 # in us
    distances = np.ones(len(parfiles)) * 1.

    PSRs = []
    for p in parfiles:
        psr = LT.fakepulsar(p, obstimes, toaerr)
        PSRs.append(psr)
     
    fgw = np.array(args.fgw)*1e-9
    mc = 10**np.array(args.lmc)

    ePSRs = prime_pulsars(PSRs, distances,
                          args.signal,
                          args.ncgw, fgw, mc, zeta = args.zeta,
                          Filter=False, psrTerm = args.psrTerm)

    with open(args.outdir+'psrs.pkl', 'wb') as psrpickle:
        pickle.dump(ePSRs, psrpickle)
    psrpickle.close()
    print('Simulation done')

    return 0
    

if __name__ == '__main__':
    main()

    

