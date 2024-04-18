# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 09:39:07 2024

@author: kgrun
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
import libstempo as T
import libstempo.plot as LP, libstempo.toasim as LT
from libstempo import spharmORFbasis as anis
import glob
import math
import json
import scipy.interpolate as interp
from enterprise.signals import gp_signals
from enterprise_extensions import model_utils, blocks
import dynesty
from enterprise.signals import signal_base
from enterprise.pulsar import Pulsar
from enterprise_extensions.frequentist import optimal_statistic as opt_stat
import enterprise.signals.parameter as parameter
from enterprise.signals import white_signals

import corner
from enterprise_extensions import sampler as sp
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc



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

def extrap1d(interpolator):
    """
    Function to extend an interpolation function to an
    extrapolation function.

    :param interpolator: scipy interp1d object

    :returns ufunclike: extension of function to extrapolation
    """

    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]  # +(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]  # +(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(list(map(pointwise, np.array(xs))))

    return ufunclike

def createGWB(
    psr,
    Amp,
    gam,
    noCorr=False,
    seed=None,
    turnover=False,
    clm=[np.sqrt(4.0 * np.pi)],
    lmax=0,
    f0=1e-9,
    beta=1,
    power=1,
    userSpec=None,
    npts=600,
    howml=1,
):
    """
    Function to create GW-induced residuals from a stochastic GWB as defined
    in Chamberlin, Creighton, Demorest, et al. (2014).

    :param psr: pulsar object for single pulsar
    :param Amp: Amplitude of red noise in GW units
    :param gam: Red noise power law spectral index
    :param noCorr: Add red noise with no spatial correlations
    :param seed: Random number seed
    :param turnover: Produce spectrum with turnover at frequency f0
    :param clm: coefficients of spherical harmonic decomposition of GW power
    :param lmax: maximum multipole of GW power decomposition
    :param f0: Frequency of spectrum turnover
    :param beta: Spectral index of power spectram for f << f0
    :param power: Fudge factor for flatness of spectrum turnover
    :param userSpec: User-supplied characteristic strain spectrum
                     (first column is freqs, second is spectrum)
    :param npts: Number of points used in interpolation
    :param howml: Lowest frequency is 1/(howml * T)

    :returns: list of residuals for each pulsar
    """

    if seed is not None:
        np.random.seed(seed)

    # number of pulsars
    Npulsars = len(psr)

    # gw start and end times for entire data set
    start = np.min([p.toas().min() * 86400 for p in psr]) - 86400
    stop = np.max([p.toas().max() * 86400 for p in psr]) + 86400

    # duration of the signal
    dur = stop - start

    # get maximum number of points
    if npts is None:
        # default to cadence of 2 weeks
        npts = dur / (86400 * 14)

    # make a vector of evenly sampled data points
    ut = np.linspace(start, stop, npts)

    # time resolution in days
    dt = dur / npts

    # compute the overlap reduction function
    if noCorr:
        ORF = np.diag(np.ones(Npulsars) * 2)
    else:
        psrlocs = np.zeros((Npulsars, 2))

        for ii in range(Npulsars):
            if "RAJ" and "DECJ" in psr[ii].pars():
                psrlocs[ii] = np.double(psr[ii]["RAJ"].val), np.double(psr[ii]["DECJ"].val)
            elif "ELONG" and "ELAT" in psr[ii].pars():
                fac = 180.0 / np.pi
                # check for B name
                if "B" in psr[ii].name:
                    epoch = "1950"
                else:
                    epoch = "2000"
                coords = ephem.Equatorial(
                    ephem.Ecliptic(str(psr[ii]["ELONG"].val * fac), str(psr[ii]["ELAT"].val * fac)), epoch=epoch
                )
                psrlocs[ii] = float(repr(coords.ra)), float(repr(coords.dec))

        psrlocs[:, 1] = np.pi / 2.0 - psrlocs[:, 1]
        anisbasis = np.array(anis.CorrBasis(psrlocs, lmax))
        ORF = sum(clm[kk] * anisbasis[kk] for kk in range(len(anisbasis)))
        ORF *= 2.0

    # Define frequencies spanning from DC to Nyquist.
    # This is a vector spanning these frequencies in increments of 1/(dur*howml).
    f = np.arange(1 / dur, 1 / (2 * dt), 1 / (dur * howml))
    f[0] = f[1]  # avoid divide by 0 warning
    Nf = len(f)

    # Use Cholesky transform to take 'square root' of ORF
    M = np.linalg.cholesky(ORF)

    # Create random frequency series from zero mean, unit variance, Gaussian distributions
    w = np.zeros((Npulsars, Nf), complex)
    for ll in range(Npulsars):
        w[ll, :] = np.random.randn(Nf) + 1j * np.random.randn(Nf)

    # strain amplitude
    if userSpec is None:

        f1yr = 1 / 3.16e7
        alpha = -0.5 * (gam - 3)
        hcf = Amp * (f / f1yr) ** (alpha)
        if turnover:
            si = alpha - beta
            hcf /= (1 + (f / f0) ** (power * si)) ** (1 / power)

    elif userSpec is not None:

        freqs = userSpec[:, 0]
        if len(userSpec[:, 0]) != len(freqs):
            raise ValueError("Number of supplied spectral points does not match number of frequencies!")
        else:
            fspec_in = interp.interp1d(np.log10(freqs), userSpec[:, 1], kind="linear")
            fspec_ex = extrap1d(fspec_in)
            hcf = 10.0 ** np.log10(userSpec[:, 1])
            #hcf = 10.0 ** np.log10(fspec_ex(f))
            print("userspec: ", np.log10(userSpec[:, 1]))
            print(hcf)
    plt.plot(hcf)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()
    plt.clf()

    C = 1 / 96 / np.pi**2 * hcf**2 / f**3 * dur * howml
    
    plt.plot(C)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()
    plt.clf()

    # inject residuals in the frequency domain
    Res_f = np.dot(M, w)
    for ll in range(Npulsars):
        Res_f[ll] = Res_f[ll] * C ** (0.5)  # rescale by frequency dependent factor
        Res_f[ll, 0] = 0  # set DC bin to zero to avoid infinities
        Res_f[ll, -1] = 0  # set Nyquist bin to zero also

    # Now fill in bins after Nyquist (for fft data packing) and take inverse FT
    Res_f2 = np.zeros((Npulsars, 2 * Nf - 2), complex)
    Res_t = np.zeros((Npulsars, 2 * Nf - 2))
    Res_f2[:, 0:Nf] = Res_f[:, 0:Nf]
    Res_f2[:, Nf : (2 * Nf - 2)] = np.conj(Res_f[:, (Nf - 2) : 0 : -1])
    Res_t = np.real(np.fft.ifft(Res_f2) / dt)

    # shorten data and interpolate onto TOAs
    Res = np.zeros((Npulsars, npts))
    res_gw = []
    for ll in range(Npulsars):
        Res[ll, :] = Res_t[ll, 10 : (npts + 10)]
        f = interp.interp1d(ut, Res[ll, :], kind="linear")
        res_gw.append(f(psr[ll].toas() * 86400))

    # return res_gw
    ct = 0
    for p in psr:
        p.stoas[:] += res_gw[ct] / 86400.0
        ct += 1

        
def createFreq(
    psr,
    seed=None,
    npts=600,
    howml=1
):
    """
    Function to create GW-induced residuals from a stochastic GWB as defined
    in Chamberlin, Creighton, Demorest, et al. (2014).

    :param psr: pulsar object for single pulsar
    :param Amp: Amplitude of red noise in GW units
    :param gam: Red noise power law spectral index
    :param noCorr: Add red noise with no spatial correlations
    :param seed: Random number seed
    :param turnover: Produce spectrum with turnover at frequency f0
    :param clm: coefficients of spherical harmonic decomposition of GW power
    :param lmax: maximum multipole of GW power decomposition
    :param f0: Frequency of spectrum turnover
    :param beta: Spectral index of power spectram for f << f0
    :param power: Fudge factor for flatness of spectrum turnover
    :param userSpec: User-supplied characteristic strain spectrum
                     (first column is freqs, second is spectrum)
    :param npts: Number of points used in interpolation
    :param howml: Lowest frequency is 1/(howml * T)

    :returns: list of residuals for each pulsar
    """

    if seed is not None:
        np.random.seed(seed)

    # number of pulsars
    Npulsars = len(psr)

    # gw start and end times for entire data set
    start = np.min([p.toas().min() * 86400 for p in psr]) - 86400
    stop = np.max([p.toas().max() * 86400 for p in psr]) + 86400

    # duration of the signal
    dur = stop - start

    # get maximum number of points
    if npts is None:
        # default to cadence of 2 weeks
        npts = dur / (86400 * 14)

    # make a vector of evenly sampled data points
    ut = np.linspace(start, stop, npts)

    # time resolution in days
    dt = dur / npts

    # Define frequencies spanning from DC to Nyquist.
    # This is a vector spanning these frequencies in increments of 1/(dur*howml).
    f = np.arange(1 / dur, 1 / (2 * dt), 1 / (dur * howml))
    f[0] = f[1]  # avoid divide by 0 warning
    Nf = len(f)
    
    return f


def run_sampler(pta, outdir = ''):

    N = int(1e7)                                    # number of samples
    x0 = np.hstack(p.sample() for p in pta.params)  # initial parameter vector
    ndim = len(x0)                                  # number of dimensions
    print('x0 =', x0)

    # initial jump covariance matrix
    cov = np.diag(np.ones(ndim) * 0.01**2)
    
    #initialize the sampler object
    sampler = ptmcmc(ndim, pta.get_lnlikelihood, pta.get_lnprior, cov, outDir=outdir, resume=False)
    
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



def produce_output(outdir = '',):
    
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
        
    with open(outdir+"/maxlike.json", "w") as f:
        json.dump(outdict, f, indent=4)   
    f.close()
    
    return chain[burn:], outdict


'''
parfiles = sorted(glob.glob('epta_init/*.par'))
Npsr = len(parfiles)
print(parfiles)


psrs = []
Amp = 2e-14
gamma = 13./3.
howml = 2

for ii in range(0,Npsr):

    # years of observations>
    psr = LT.fakepulsar(parfile=parfiles[ii],
            obstimes=np.arange(53000,53000+10*365.25,28.), toaerr=0.1)

    # We now remove the computed residuals from the TOAs, obtaining (in effect) a perfect realization of the deterministic timing model. The pulsar parameters will have changed somewhat, so `make_ideal` calls `fit()` on the pulsar object.
    LT.make_ideal(psr)

    #Generate white noise
    LT.add_efac(psr,efac=1.0)

    # add to list
    psrs.append(psr)

#LT.createGWB(psrs, Amp=Amp, gam=gamma, howml=100)

#make spectrum
#define the spectrum
freq = createFreq(psrs, howml=howml)
spec = 1e-50*np.ones(len(freq))
spec[7*howml] = 1e-12*np.ones(1)
userSpec = np.asarray([freq, spec]).T

createGWB(psrs, Amp=Amp, gam=gamma, howml=howml, userSpec=userSpec)

for Psr in psrs:

    Psr.savepar("epta_sim/" + Psr.name + '.par')
    Psr.savetim("epta_sim/" + Psr.name + '.tim')
    T.purgetim("epta_sim/" + Psr.name + '.tim')
'''


parfiles = sorted(glob.glob('epta_sim/*.par'))
timfiles = sorted(glob.glob('epta_sim/*.tim'))

psrs = []
ephemeris = None
for p, t in zip(parfiles, timfiles):
    psr = Pulsar(p, t, ephem=ephemeris)
    psrs.append(psr)

pickle.dump(psrs, open('psrs.pkl', 'wb'))

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
s += blocks.common_red_noise_block(psd='spectrum', prior='log-uniform', Tspan=Tspan,
                                   components=30, name='gw_crn', orf = None)

# We set up the PTA object using the signal we defined above and the pulsars
pta = signal_base.PTA([s(p) for p in psrs]) 
    
print(pta.params)

run_sampler(pta, "./epta_sim/")

fs = (np.arange(30) + 1) / Tspan
chainname = 'chain_1'
chain = np.loadtxt("./epta_sim/" + chainname + '.txt')


chain, setpars = produce_output("./epta_sim/")

ostat = opt_stat.OptimalStatistic(psrs, pta=pta, orf='hd')

print('OS type: free spectrum')
print(50*'-' + '\n' + 50*'-')
print('Calculating OS stats (maximum LH)')

out_maxLH = np.zeros((len(fs), 4))

for ii, f in enumerate(fs):
    (10*'-', f)        
    xi, rho, sig, OS, OS_sig = ostat.compute_os(params=setpars, psd='spectrum', fgw=f)
    xi_mean, xi_err, rho_avg, sig_avg = sort_and_bin(xi, rho, sig)
    OS_output = np.array([f, OS, OS_sig, 0])
    HD_output = np.array([xi_mean, xi_err, np.array(rho_avg), np.array(sig_avg)])
    plt.errorbar(xi_mean*180/np.pi, rho_avg, xerr=xi_err*180/np.pi, marker='o', ls='', color='k', capsize=4)
    eta = np.linspace(0.01,180,100)
    HD = get_HD_curve(eta+1)
    plt.plot(eta, OS*HD, ls='--', label='Hellings-Downs', color='b')
    plt.savefig('./epta_sim/HDcurve_{}.png'.format('{:.2e}'.format(f)))
    plt.clf()
        
    #np.savetxt('./output/HD_output_WNRN_{:.2e}.txt'.format(f), HD_output.transpose(), delimiter='\t')
    out_maxLH[ii] = OS_output
    np.savetxt('./epta_sim/GWBsinglebin.txt', out_maxLH, delimiter='\t')
    print()
    print(50*'-' + '\n' + 50*'-')
    print('Calculating OS stats (noise draws)')
    N = 1000
    print('...Iterating over {} chain draws'.format(N))
    out_noisedraws = np.zeros((len(fs)*3, N))
    
    for ii, f in enumerate(fs):
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
        
np.savetxt('./epta_sim/GWBsinglebin_NM.txt', out_noisedraws, delimiter='\t') 


