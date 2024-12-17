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
                LT.add_rednoise(ltp,10**np.random.uniform(-14.,-13.),np.random.uniform(1.7,2.2), components=20)
                outstr += ' added RN'
            
            elif s=='CRN':
                LT.add_rednoise(ltp,10**-14,2.2, tspan=3*24*3600, components=10)
                outstr += ' added CRN'
 
            elif s=='GWB':
                gamma = 13./3.
                alpha = (3.-gamma)/2.
                log10_A_gw = -14.5
                flow = 1e-9
                fhigh = 1e-8
                LT.add_gwb(ltp, gwAmp=10**log10_A_gw, alpha=alpha, flow=flow, fhigh=fhigh)
                outstr += ' added GWB with log10A={}, alpha={:.2}'.format(log10_A_gw,alpha)
            
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
                
                if Filter == True:
                    if f0 - 5/(3652*24*3600) > 1/(3652*24*3600):
                        low = f0 - 5/(3652*24*3600)
                    else:
                        low = 1/(3652*24*3600)


                    if f0 + 5/(3652*24*3600) < 1/(2*14*24*3600):
                        high = f0 + 5/(3652*34*3600)
                    else:
                        high = 1/(1.5*14*24*3600)

                    sos = signal.butter(20, [low, high], btype='bandpass', output='sos', fs = 1/(14*24*3600), analog=False)
                    filtered_res = signal.sosfiltfilt(sos, ltp.residuals())

                    if f0 == 16/(3652*24*3600):
                        ltp.stoas[:] = stoas_temp + filtered_res/86400
                        
                        plt.plot(ltp.toas(), filtered_res)
                        plt.plot(ltp.toas(), ltp.residuals())
                        plt.savefig('./{}.png'.format(ltp.name))
                        plt.clf()
            
            else:
                print('unsupported signal type')
                break
            
        print(outstr)
    
        ePSRs.append(Pulsar(ltp, dist=pdistances[ii]))
        
    return ePSRs



def PTA_model(psrs, obtimes, ptamodel, mode=None):
    #find the maximum time span to set GW frequency sampling
    Tspan = model_utils.get_tspan(psrs)
    
    s = []
    model = ptamodel.split(',')
    outstring = 'PTA model with '
    
    for m in model:
        if m =='TM':
            s.append(gp_signals.TimingModel())   # First we add the timing model
            outstring += 'TM '

        elif m=='WN':
            s.append(white_signals.MeasurementNoise(efac=1.))   #Timing noise
            outstring += 'WN '
    
        elif m=='RN':
            #log10_A = parameter.Uniform(-15., -11.)
            #gamma = parameter.Uniform(1.5, 2.5)
            #pl = utils.powerlaw(log10_A=log10_A, gamma=gamma)
            #s.append(gp_signals.FourierBasisGP(spectrum=pl, Tspan=Tspan, components=20) )
            s.append(blocks.red_noise_block(components = 20))
            outstring += 'RN '
        
        elif m=='CRN_fs':
            s.append(blocks.common_red_noise_block(psd='spectrum', prior='log-uniform', name='gw', Tspan = Tspan, components=20))
            outstring += 'CRN (spectrum) '

        elif m=='CRN_pl':
            s.append(blocks.common_red_noise_block(psd='powerlaw', prior='log-uniform', Tspan = Tspan, name='gw', components = 20))
            outstring += 'CRN (powerlaw) '

        else:
            print('Unsupported PTA model component')
            break
    
    signal = s[0]
    for sl in s[1:]:
        signal+=sl
    
    print(outstring)
   

    # We set up the PTA object using the signal we defined above and the pulsars
    pta = signal_base.PTA([signal(p) for p in psrs])
    
    return pta



def run_sampler(pta, outdir = '', resume=False):

    N = int(1e6)                                    # number of samples
    x0 = np.hstack(p.sample() for p in pta.params)  # initial parameter vector
    ndim = len(x0)                                  # number of dimensions
    print('x0 =', x0)

    # initial jump covariance matrix
    cov = np.diag(np.ones(ndim) * 0.01**2)
    
    #initialize the sampler object
    sampler = ptmcmc(ndim, pta.get_lnlikelihood, pta.get_lnprior, cov, outDir=outdir, resume=resume)
    
    # additional jumps
    jp = sp.JumpProposal(pta)
    sampler.addProposalToCycle(jp.draw_from_prior, 5)
    
    sel_sig = ["rn", "red_noise", "dm_gp", "fcn", "chrom-rn", "srn", "dm_srn", "freechrom-srn", "chrom-srn",
                        "dm-expd", "freechrom-expd", "chrom-expd",
                        "dm-y", "freechrom-y", "chrom-y",
                        "gw"]
    for s in sel_sig:
        if any([s in p for p in pta.param_names]):
            #pnames = [p.name for p in pta.params if s in p.name]
            #print('Adding %s prior draws with parameters :'%s, pnames, '\n')
            print('Adding %s prior draws.'%s)
            sampler.addProposalToCycle(jp.draw_from_par_prior(s), 10)

        
    sampler.sample(x0, N, SCAMweight=40, AMweight=25, DEweight=55) # these weights relate to frequency of jumps

    # write a list of the parameters to a text file
    # and a list of the parameter groupings used
    filename = outdir + '/params.txt'
    np.savetxt(filename,list(map(str, pta.param_names)), fmt='%s')
    
    return None




def produce_output(outdir = '', cornerplot=False):
    
    chainname = 'chain_1'
    chain = np.loadtxt(outdir + chainname + '.txt')
    
    params = []
    with open(outdir + 'params.txt', 'r') as file:
        for line in file.readlines():
            params.append(line.strip())
    


    burn = int(0.3*chain.shape[0])   # burn = beginning of chain, not helpful!
    outdict={}
    
    for i, p in enumerate(params):    
        # sigma
        lisst = corner.quantile(chain[burn:,i], [0.16, 0.5, 0.84])
        # replace mean with maximum likelihood
        n, bins, patches = plt.hist(chain[burn:,i], bins=100)
        plt.clf()
        max_lh_pos = np.argmax(n)
        max_lh = bins[max_lh_pos]
        lisst[1] = max_lh
        
        outdict[p]=max_lh
        
    if cornerplot==True:
        corner.corner(chain[burn:,:-4], labels=params[:],
                      bins =30,
                      plot_datapoints=False, plot_density=True, 
                      plot_contours=False,fill_contours=False,
                      show_titles = True, use_math_text=True, quantiles=[0.16, 0.5, 0.84], verbose=True)
        plt.savefig('{}/cornerplot.png'.format(outdir), bbox_inches='tight')
        plt.clf()
        
    with open(outdir+"/maxlike.json", "w") as f:
        json.dump(outdict, f, indent=4)   
    f.close()
    
    return chain[burn:]
    


# Before plotting, we need to bin the cross-correlations

def weightedavg(rho, sig):
    weights, avg = 0., 0.
    for r,s in zip(rho,sig):
        weights += 1./(s*s)
        avg += r/(s*s)
        
    return avg/weights, np.sqrt(1./weights)



def sort_and_bin(xi, rho, sig):
    # sort the cross-correlations by xi
    idx = np.argsort(xi)

    xi_sorted = xi[idx]
    rho_sorted = rho[idx]
    sig_sorted = sig[idx]

    # bin the cross-correlations so that there are the same number of pairs per bin
    npairs = 10

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
    
    return xi_mean, xi_err, rho_avg, sig_avg
    
    
    
    
def get_HD_curve(zeta):
    coszeta = np.cos(zeta*np.pi/180.)
    xip = (1.-coszeta) / 2.
    HD = 3.*( 1./3. + xip * ( np.log(xip) -1./6.) )
    
    return HD/2
    



def do_analysis(PSRs, obstimes, signal, ptamodel, Ncgw = 1, fgw=[5e-8], mc=[10**9], zeta=1., outdir = '', cornerplot=False, Filter=True, psrTerm = True, resume=False):
    #distances = np.random.uniform(0.8,1.2, len(PSRs))
    distances = np.ones(len(PSRs)) * 1.
    
    if isinstance(PSRs, list):
        ePSRs = prime_pulsars(PSRs, distances, signal, Ncgw, fgw, mc, zeta = zeta, Filter=Filter, psrTerm=psrTerm)
    else:
        with open(PSRs, 'rb') as psrpickle:
            ePSRs = pickle.load(psrpickle)
        print('loaded pulsars from pickle')
        psrpickle.close()

    pta = PTA_model(ePSRs, obstimes, ptamodel)

    run_sampler(pta, '{}/'.format(outdir), resume = resume)
    chain = produce_output('{}/'.format(outdir), cornerplot=cornerplot)
        
    with open('{}/maxlike.json'.format(outdir), 'r') as f:
        ml_params = json.load(f)
    f.close()
        
        
    setpars = (pta.map_params(list(ml_params.values())))
    setpars.update(ml_params)
        
    return ePSRs, pta, setpars, chain





def set_up_global_options():
    parser = argparse.ArgumentParser(description='Generate sensitivity curves for PTA datasets.')
    parser.add_argument('--par', type=str, default=None, help='Path to the parfiles.')
    parser.add_argument('--outdir', type=str, default=None, help='Directory for the MCMC output')
    parser.add_argument('--result', type=str, default=None, help='Directory and filename for the OS stats output')

    parser.add_argument('--signal', type=str, default='CGW', help='Which signals to inject to the pulsar objects')
    parser.add_argument('--ptamodel', type=str, default='CGW', help='Which signals to include in the PTA model')
    parser.add_argument('--lmc', type=float, nargs='*', default=9.5)
    parser.add_argument('--fgw', type=float, nargs='*', default=22.3)
    parser.add_argument('--ncgw', type=int, default=1)

    parser.add_argument('--psrTerm', action="store_true", default=False, help='Simulate CGW with or without pulsar term')
    parser.add_argument('--zeta', type=float, default=1.0)
    parser.add_argument('--cornerplot', action="store_true", default=False, help='Save cornerplot from the MCMC run')
    parser.add_argument('--N', type=int, default=None, help='number of MCMC draws')
    parser.add_argument('--noisemarginalised', action="store_true", default=False, help='draw samples from MCMC to get distribution of OS values')
    parser.add_argument('--psrpickle', type=str, default=None, help='psrpickle')
    parser.add_argument('--OStype', type=str, default='spectrum', help='type for OS')
    parser.add_argument('--resume', action="store_true", default=False, help='Resume MCMC run')

    return parser.parse_args()






    
    
def main():
    args = set_up_global_options()
    
    ### load pulsars ###
    obstimes = np.arange(50000,53652,14)
    if args.psrpickle is not None:
        PSRs = args.psrpickle

    else:
        parfiles = glob.glob(args.par+'/*.par')

        toaerr = 0.1 # in us

        PSRs = []
        for p in parfiles:
            psr = LT.fakepulsar(p, obstimes, toaerr)
            PSRs.append(psr)


    ePSRs, pta, setpars, chain = do_analysis(PSRs, obstimes, args.signal, args.ptamodel, fgw=np.array(args.fgw)*1e-9, mc=10**np.array(args.lmc), zeta=args.zeta, outdir=args.outdir, Ncgw=args.ncgw,  cornerplot=args.cornerplot, psrTerm=args.psrTerm, Filter=False, resume = args.resume)

    tspan = model_utils.get_tspan(ePSRs)
    tspan_yrs = tspan/(3600*24*365.25)
    print('This PTA has a total time span of {} years \n'.format(tspan_yrs))
    print('Initiating OS stats')
    ostat = opt_stat.OptimalStatistic(ePSRs, pta=pta, orf='hd')

    if args.OStype == 'spectrum':

        test_frequencies = np.linspace(1/tspan, 20/tspan, 20) 
        #test_frequencies = np.linspace(1/(3625*86400), 20/(3625*86400), 20)

        print('OS type: free spectrum')
        print(50*'-' + '\n' + 50*'-')
        print('Calculating OS stats (maximum LH)')

        out_maxLH = np.zeros((len(test_frequencies), 4))

        for ii, f in enumerate(test_frequencies):
            (10*'-', f)
            xi, rho, sig, OS, OS_sig = ostat.compute_os(params=setpars, psd='spectrum', fgw=f)
            xi_mean, xi_err, rho_avg, sig_avg = sort_and_bin(xi, rho, sig)

            OS_output = np.array([f, OS, OS_sig, args.zeta])
            HD_output = np.array([xi_mean, xi_err, np.array(rho_avg), np.array(sig_avg)])

            plt.errorbar(xi_mean*180/np.pi, rho_avg, xerr=xi_err*180/np.pi, marker='o', ls='', color='k', capsize=4)
            eta = np.linspace(0.01,180,100)
            HD = get_HD_curve(eta+1)
            plt.plot(eta, OS*HD, ls='--', label='Hellings-Downs', color='b')
            plt.savefig('{}/HDcurve_{}.png'.format(args.outdir, '{:.2e}'.format(f)))
            plt.clf()

            #np.savetxt('./output/HD_output_WNRN_{:.2e}.txt'.format(f), HD_output.transpose(), delimiter='\t')
            out_maxLH[ii] = OS_output

        np.savetxt(args.result + '.txt', out_maxLH, delimiter='\t')
        '''
        print()
        print(50*'-' + '\n' + 50*'-')
        print('Calculating OS stats (noise draws)')
        N = args.N
        print('...Iterating over {} chain draws'.format(N))
        out_noisedraws = np.zeros((len(test_frequencies)*3, N))
        
        for ii, f in enumerate(test_frequencies):
            print('-- {}'.format(f), flush=True)
            opt_temp = []
            sig_temp = []
            for n in range(N):
                idx = np.random.randint(0, chain.shape[0])

                print('...on sample {}'.format(n), end='\r', flush=True)

                for jj, p in enumerate(pta.param_names):
                    setpars[p] = chain[idx][jj]
                    setpars['gw_log10_rho'] = chain[idx][-24:-4]

                # compute OS stat for the drawn parameters
                _, _, _, OS_temp, OS_sig_temp = ostat.compute_os(params=setpars, psd='spectrum', fgw=f)

                opt_temp.append(OS_temp)
                sig_temp.append(OS_sig_temp)

            opt = np.array(opt_temp)
            sig = np.array(sig_temp)
            sn = opt/sig

            out_noisedraws[int(3*ii)] = opt
            out_noisedraws[int(3*ii+1)] = sig
            out_noisedraws[int(3*ii+2)] = sn

        np.savetxt(args.result + '_NM.txt', out_noisedraws, delimiter='\t')
        '''

    #------------------------------------------------------------
   
    elif args.OStype == 'powerlaw':
        print('OS type: powerlaw')
        print(50*'-' + '\n' + 50*'-')
        print('Calculating OS stats (maximum LH)')

        xi, rho, sig, OS, OS_sig = ostat.compute_os(params=setpars, psd='powerlaw')
        out_maxLH = np.array([OS, OS_sig, OS/OS_sig])
        print(out_maxLH)


        print(50*'-' + '\n' + 50*'-')
        print('Calculating OS stats (noise draws)')
        N = args.N
        print('...Iterating over {} chain draws'.format(N))
        out_noisedraws = np.zeros((3, N))

        opt_temp = []
        sig_temp = []
        for n in range(N):
            idx = np.random.randint(0, chain.shape[0])

            print('...on sample {}'.format(n), end='\r', flush=True)

            for jj, p in enumerate(pta.param_names):
                setpars[p] = chain[idx][jj]
                #setpars['gw_log10_rho'] = chain[idx][:-4]

            # compute OS stat for the drawn parameters
            _, _, _, OS_temp, OS_sig_temp = ostat.compute_os(params=setpars, psd='powerlaw')

            opt_temp.append(OS_temp)
            sig_temp.append(OS_sig_temp)

        opt = np.array(opt_temp)
        sig = np.array(sig_temp)
        sn = opt/sig

        out_noisedraws[0] = opt
        out_noisedraws[1] = sig
        out_noisedraws[2] = sn

        np.savetxt(args.result+'_NM.txt', out_noisedraws.transpose(), delimiter='\t')
 
    return 0
    

if __name__ == '__main__':
    main()

    

