# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 16:48:15 2024

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

import ephem
import acor

import libstempo.toasim as LT
import libstempo.plot as LP

from enterprise import constants as const
from enterprise.signals import signal_base, gp_signals, white_signals, utils, gp_priors
from enterprise.pulsar import Pulsar
import enterprise.signals.parameter as parameter

from enterprise_extensions import model_utils, blocks, deterministic
from enterprise_extensions.frequentist import optimal_statistic as opt_stat
from enterprise_extensions import sampler

import corner
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc
from chainconsumer import ChainConsumer


@signal_base.function
def cw_delay(toas, pos, pdist,
             cos_gwtheta=0, gwphi=0, cos_inc=0,
             log10_mc=9, log10_fgw=-8, log10_dist=None, log10_h=None,
             phase0=0, psi=0,
             psrTerm=False, p_dist=None, p_phase=None,
             evolve=False, phase_approx=False, check=False,
             tref=0):
    """
    Function to create GW incuced residuals from a SMBMB as
    defined in Ellis et. al 2012,2013.

    :param toas:
        Pular toas in seconds
    :param pos:
        Unit vector from the Earth to the pulsar
    :param pdist:
        Pulsar distance (mean and uncertainty) [kpc]
    :param cos_gwtheta:
        Cosine of Polar angle of GW source in celestial coords [radians]
    :param gwphi:
        Azimuthal angle of GW source in celestial coords [radians]
    :param cos_inc:
        cosine of Inclination of GW source [radians]
    :param log10_mc:
        log10 of Chirp mass of SMBMB [solar masses]
    :param log10_fgw:
        log10 of Frequency of GW (twice the orbital frequency) [Hz]
    :param log10_dist:
        log10 of Luminosity distance to SMBMB [Mpc],
        used to compute strain, if not None
    :param log10_h:
        log10 of GW strain,
        used to compute distance, if not None
    :param phase0:
        Initial GW phase of source [radians]
    :param psi:
        Polarization angle of GW source [radians]
    :param psrTerm:
        Option to include pulsar term [boolean]
    :param p_dist:
        Pulsar distance parameter
    :param p_phase:
        Use pulsar phase to determine distance [radian]
    :param evolve:
        Option to include/exclude full evolution [boolean]
    :param phase_approx:
        Option to include/exclude phase evolution across observation time
        [boolean]
    :param check:
        Check if frequency evolves significantly over obs. time [boolean]
    :param tref:
        Reference time for phase and frequency [s]

    :return: Vector of induced residuals

    """

    # convert units to time
    mc = 10**log10_mc * const.Tsun
    fgw = 10**log10_fgw
    gwtheta = np.arccos(cos_gwtheta)
    inc = np.arccos(cos_inc)
    
    p_dist = p_dist*const.kpc/const.c

    if log10_h is None and log10_dist is None:
        raise ValueError("one of log10_dist or log10_h must be non-None")
    elif log10_h is not None and log10_dist is not None:
        raise ValueError("only one of log10_dist or log10_h can be non-None")
    elif log10_h is None:
        dist = 10**log10_dist * const.Mpc / const.c
    else:
        dist = 2 * mc**(5/3) * (np.pi*fgw)**(2/3) / 10**log10_h

    if check:
        # check that frequency is not evolving significantly over obs. time
        fstart = fgw * (1 - 256/5 * mc**(5/3) * fgw**(8/3) * toas[0])**(-3/8)
        fend = fgw * (1 - 256/5 * mc**(5/3) * fgw**(8/3) * toas[-1])**(-3/8)
        df = fend - fstart

        # observation time
        Tobs = toas.max()-toas.min()
        fbin = 1/Tobs

        if np.abs(df) > fbin:
            print('WARNING: Frequency is evolving over more than one '
                  'frequency bin.')
            print('f0 = {0}, f1 = {1}, df = {2}, fbin = {3}'.format(fstart, fend, df, fbin))
            return np.ones(len(toas)) * np.nan

    # get antenna pattern funcs and cosMu
    # write function to get pos from theta,phi
    fplus, fcross, cosMu = utils.create_gw_antenna_pattern(pos, gwtheta, gwphi)

    # get pulsar time
    toas -= tref
    if p_dist > 0:
        tp = toas-p_dist*(1-cosMu)
    else:
        tp = toas

    # orbital frequency
    w0 = np.pi * fgw
    phase0 /= 2  # convert GW to orbital phase
    # omegadot = 96/5 * mc**(5/3) * w0**(11/3) # Not currently used in code

    # evolution
    if evolve:
        # calculate time dependent frequency at earth and pulsar
        #omega = w0 * (1 - 256/5 * mc**(5/3) * w0**(8/3) * toas)**(-3/8)
        #omega_p = w0 * (1 - 256/5 * mc**(5/3) * w0**(8/3) * tp)**(-3/8)

        omega = w0 * np.sign( (1 - 256/5 * mc**(5/3) * w0**(8/3) * toas)) * (np.abs((1 - 256/5 * mc**(5/3) * w0**(8/3) * toas))) ** (-3/8)
        omega_p = w0 * np.sign( (1 - 256/5 * mc**(5/3) * w0**(8/3) * toas)) * (np.abs((1 - 256/5 * mc**(5/3) * w0**(8/3) * toas))) ** (-3/8)


        if p_dist > 0:
            omega_p0 = w0 * (1 + 256/5
                             * mc**(5/3) * w0**(8/3) * p_dist*(1-cosMu))**(-3/8)
        else:
            omega_p0 = w0

        # calculate time dependent phase
        #phase = phase0 + 1/32/mc**(5/3) * (w0**(-5/3) - omega**(-5/3))
        phase = phase0 + 1/32/mc**(5/3) * (w0**(-5/3) - np.sign(omega) * np.abs(omega)**(-5/3))

        if p_phase is None:
            phase_p = phase0 + 1/32/mc**(5/3) * (w0**(-5/3) - np.sign(omega_p) * np.abs(omega_p)**(-5/3))
        else:
            phase_p = (phase0 + p_phase
                       + 1/32*mc**(-5/3) * (omega_p0**(-5/3) - np.sign(omega_p)*np.abs(omega_p)**(-5/3)))

    elif phase_approx:
        # monochromatic
        omega = w0
        if p_dist > 0:
            omega_p = w0 * (1 + 256/5
                            * mc**(5/3) * w0**(8/3) * p_dist*(1-cosMu))**(-3/8)
        else:
            omega_p = w0

        # phases
        phase = phase0 + omega * toas
        if p_phase is not None:
            phase_p = phase0 + p_phase + omega_p*toas
        else:
            phase_p = (phase0 + omega_p*toas
                       + 1/32/mc**(5/3) * (w0**(-5/3) - omega_p**(-5/3)))

    # no evolution
    else:
        # monochromatic
        omega = np.pi*fgw
        omega_p = omega

        # phases
        phase = phase0 + omega * toas
        phase_p = phase0 + omega * tp

    # define time dependent coefficients
    At = -0.5*np.sin(2*phase)*(3+np.cos(2*inc))
    Bt = 2*np.cos(2*phase)*np.cos(inc)
    At_p = -0.5*np.sin(2*phase_p)*(3+np.cos(2*inc))
    Bt_p = 2*np.cos(2*phase_p)*np.cos(inc)

    # now define time dependent amplitudes
    alpha = mc**(5./3.)/(dist*omega**(1./3.))
    alpha_p = mc**(5./3.)/(dist*omega_p**(1./3.))

    # define rplus and rcross
    rplus = alpha*(-At*np.cos(2*psi)+Bt*np.sin(2*psi))
    rcross = alpha*(At*np.sin(2*psi)+Bt*np.cos(2*psi))
    rplus_p = alpha_p*(-At_p*np.cos(2*psi)+Bt_p*np.sin(2*psi))
    rcross_p = alpha_p*(At_p*np.sin(2*psi)+Bt_p*np.cos(2*psi))

    # residuals
    if psrTerm:
        res = fplus*(rplus_p-rplus)+fcross*(rcross_p-rcross)
    else:
        res = -fplus*rplus - fcross*rcross

    return res





def set_up_global_options():
    parser = argparse.ArgumentParser(description='Generate sensitivity curves for PTA datasets.')
    parser.add_argument('--basedir', type=str, default=None, help='Path to the parfiles.')
    parser.add_argument('--outdir', type=str, default=None, help='Path to the parfiles.')
    parser.add_argument('--noisefile', type=str, default=None, help='json file with noise parameters')

    
    parser.add_argument('--ptamodel', type=str, default='CGW', help='Which signals to include in the PTA model')
    parser.add_argument('--lmc', type=float, nargs='*', default=9.5)
    parser.add_argument('--fgw', type=float, nargs='*', default=22.3)
    parser.add_argument('--ncgw', type=int, default=1)

    parser.add_argument('--psrTerm', action="store_true", default=False, help='Simulate CGW with or without pulsar term')
    parser.add_argument('--pd', type=float, default=1.0)

    parser.add_argument('--run', action="store_true", default=False, help='run the MCMC')
    parser.add_argument('--analysis', action="store_true", default=False, help='Save cornerplot from the MCMC run')
    parser.add_argument('--use_distance', action="store_true", default=False, help='Use the randomly drawn distance')
    parser.add_argument('--sample_pdist', action="store_true", default=False, help='Use the randomly drawn distance')

    parser.add_argument('--Nsample', type=int, default=1e6, help='number of MCMC draws')
    return parser.parse_args()


def PTA_model(psrs, ptamodel, psrTerm=True, pdist=1., sample_pdist=False):
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
            efac = parameter.Constant()
            equad = parameter.Constant()

            ef = white_signals.MeasurementNoise(efac=efac)
            eq = white_signals.TNEquadNoise(log10_tnequad=equad)
            s.append(ef)   
            s.append(eq)
            outstring += 'WN '
    
        elif m=='RN':
            log10_A = parameter.Uniform(-15., -12.)
            gamma = parameter.Uniform(1.5, 2.5)
            pl = utils.powerlaw(log10_A=log10_A, gamma=gamma)
            s.append(gp_signals.FourierBasisGP(spectrum=pl, Tspan=Tspan, components=10) )
            outstring += 'RN '
        
        elif m=='CRN':
            s.append(blocks.common_red_noise_block(psd='spectrum', prior='log-uniform', name='gw', Tspan = Tspan, components=20))
            outstring += 'CRN (spectrum) '
        
        elif m=='CGW':

            cos_gwtheta = parameter.Uniform(-1, 1)('cos_gwtheta')   # position of source
            gwphi = parameter.Uniform(0, 2*np.pi)('gwphi')          # position of source
            log10_mc = parameter.Uniform(8., 10.)('log10_Mc')      # chirp mass
            log10_h = parameter.Uniform(-16, -11)('log10_h')        # strain amplitude
            log10_fgw = parameter.Uniform(-8.5, -6.5)('log10_fgw')    # gw frequency
            phase0 = parameter.Uniform(0, 2*np.pi)('phase0')        # gw phase
            psi = parameter.Uniform(0, np.pi)('psi')                # gw polarization 
            cos_inc = parameter.Uniform(-1, 1)('cos_inc')           # inclination of binary with respect to Earth

            if psrTerm == True:
                if pdist == None:
                    if sample_pdist == True:
                       p_dist = parameter.Normal(0,1)
                       #p_dist = parameter.Uniform(0,1)
                    else:
                       p_dist = 0
                    
                    p_phase = parameter.Uniform(0, np.pi)

                    cw_wf = deterministic.cw_delay(cos_gwtheta=cos_gwtheta, gwphi=gwphi, log10_mc=log10_mc,
                                                   log10_h=log10_h, log10_fgw=log10_fgw, phase0=phase0,
                                                   psi=psi, cos_inc=cos_inc,
                                                   evolve=True, psrTerm=True, phase_approx=False,
                                                   p_dist=p_dist, p_phase=p_phase)

                else:
                    p_phase = parameter.Uniform(0, np.pi)
                    p_dist = parameter.Uniform(0.7*pdist, 1.3*pdist)
                    #p_dist = pdist
                    #p_phase = None

                    cw_wf = cw_delay(cos_gwtheta=cos_gwtheta, gwphi=gwphi, log10_mc=log10_mc,
                                     log10_h=log10_h, log10_fgw=log10_fgw, phase0=phase0,
                                     psi=psi, cos_inc=cos_inc,
                                     evolve=False, psrTerm=True, phase_approx=True,
                                     p_dist=p_dist, p_phase=p_phase)


            else:
                p_phase = None
                p_dist = 0
            
                cw_wf = cw_delay(cos_gwtheta=cos_gwtheta, gwphi=gwphi, log10_mc=log10_mc, 
                                 log10_h=log10_h, log10_fgw=log10_fgw, phase0=phase0, 
                                 psi=psi, cos_inc=cos_inc,
                                 evolve = False, psrTerm=False, phase_approx=False,
                                 p_dist=p_dist, p_phase=p_phase)
            

            CGW = deterministic.CWSignal(cw_wf, ecc=False, psrTerm=psrTerm)
            #CGW = deterministic.cw_block_circ(psrTerm=True)
            s.append(CGW)
            outstring += 'CGW '
            
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




def produce_output(outdir = '' , plot = 'chainconsumer'):
    chain = np.loadtxt(outdir + '/chain_1.txt')
    
    params = []
    with open(outdir + '/pars.txt', 'r') as file:
        for line in file.readlines():
            params.append(line.strip())
    


    burn = int(0.3*chain.shape[0])   # burn = beginning of chain, not helpful!
    outdict={}
        
    if plot == 'corner':
        for i, p in enumerate(params):    
            # sigma
            parameter_estimate_list = corner.quantile(chain[burn:,i], [0.16, 0.5, 0.84])
            # replace mean with maximum likelihood
            n, bins, patches = plt.hist(chain[burn:,i], bins=100)
            plt.clf()
            max_lh_pos = np.argmax(n)
            max_lh = bins[max_lh_pos]
            parameter_estimate_list[1] = max_lh
        
            outdict[p]=parameter_estimate_list.tolist()

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
    
    elif plot == 'chainconsumer':
        cc = ChainConsumer()
        clr = "#5F9EA0" #cadetblue
        #clr = "#6199AE" #steel blue

        cc.add_chain(chain[:,:-4], parameters = params, name='cornerplot')

        cc.configure(max_ticks = 3, colors=clr,
                     tick_font_size=14, label_font_size=8, spacing=1.0,diagonal_tick_labels=True,
                     contour_labels='confidence', contour_label_font_size=14, #shade_gradient=[3.0], 
                     sigmas=[1,2,3], shade_alpha=0.6, linewidths=1.5,
                     summary=False, sigma2d=True)

        
        fig = cc.plotter.plot(legend=False, filename='{}/cornerplot.png'.format(outdir))
        fig.savefig('{}/cornerplot.png'.format(outdir), bbox_inches='tight')

        summary = cc.analysis.get_summary()
        with open(outdir+"/maxlike.json", "w") as f:
            json.dump(summary, f, indent=4)
        f.close()

    return chain[burn:]




def main():
    args = set_up_global_options()
    
    outD = args.outdir
    
    with open(args.basedir + '/psrs.pkl', 'rb') as psrpickle:
        ePSRs = pickle.load(psrpickle)
        print('loaded pulsars from pickle')
    psrpickle.close()
   
    print('...loading noise parameters')
    # load all noise parameters into a dictionary
    with open(args.noisefile, 'r') as nf:
        noisedict = json.load(nf)
    nf.close()

    if os.path.isfile(args.basedir + '/distances.json') and args.use_distance == True:
        pdistances = json.load(open(args.basedir + '/distances.json', 'r'))
        for p in ePSRs:
            p._pdist = (pdistances[p.name], 0.2)
            print(p.pdist)
        PTA = PTA_model(ePSRs, args.ptamodel, psrTerm=args.psrTerm, pdist=None, sample_pdist = args.sample_pdist)
    else:
        PTA = PTA_model(ePSRs, args.ptamodel, psrTerm=args.psrTerm, pdist=args.pd)
    
    PTA.set_default_params(noisedict)

    #  -- helpful PTA summary output
    with open(args.outdir+'/summary.txt', 'w') as sumfile:
        sumfile.write(PTA.summary())
    sumfile.close()

    print(PTA.params) 
    x0 = np.hstack([p.sample() for p in PTA.params])
    iteration = 0
    
    #print('... looking for a suitable initial sample, iteration {}'.format(iteration), end = '\r', flush=True)
    #while PTA.get_lnlikelihood(x0) < 0. or PTA.get_lnprior(x0) == float(-np.inf):
    #    iteration += 1
    #    x0 = np.hstack([p.sample() for p in PTA.params])
    #    print('\r ... looking for a suitable initial sample, iteration {}'.format(iteration), end = '\r',flush=True)
    #print(' \n ... found a sample' , flush = True)

    if args.run == True:
        smpl = sampler.setup_sampler(PTA, outdir=outD)
        smpl.sample(x0, args.Nsample, SCAMweight=60, AMweight=0, DEweight=30)#, burn=int(0.3*args.Nsample))
    
    
    if args.analysis == True:
        _ = produce_output(outdir = outD)
        #_ = produce_output(outdir = outD, plot = 'corner')
        return None
    
    else:
        return None
    

if __name__ == '__main__':
    main()
        
