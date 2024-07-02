# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 10:45:55 2024

@author: kgrun
"""


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



'''
parfiles = sorted(glob.glob('/u/kgrunthal/HD/epta_sim/par/*.par'))
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
spec[7*howml] = 5e-12*np.ones(1)
userSpec = np.asarray([freq, spec]).T

createGWB(psrs, Amp=Amp, gam=gamma, howml=howml, userSpec=userSpec)

for Psr in psrs:

    Psr.savepar("/u/kgrunthal/HD/epta_sim/par/" + Psr.name + '.par')
    Psr.savetim("/u/kgrunthal/HD/epta_sim/tim/" + Psr.name + '.tim')
    T.purgetim("/u/kgrunthal/HD/epta_sim/tim/" + Psr.name + '.tim')

'''

parfiles = sorted(glob.glob('/u/kgrunthal/HD/epta_sim/par/*.par'))
timfiles = sorted(glob.glob('/u/kgrunthal/HD/epta_sim/tim/*.tim'))

psrs = []
ephemeris = None
for p, t in zip(parfiles, timfiles):
    psr = Pulsar(p, t, ephem=ephemeris)
    psrs.append(psr)
    
    
for j in psrs:
    plt.plot(j.residuals)

plt.savefig("/u/kgrunthal/HD/epta_sim/test_1/residuals.png")
plt.clf()

    
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





def run_sampler(pta, outdir = ''):

    N = int(1e5)                                    # number of samples
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

#print(pta.param_names)

#run_sampler(pta, "/u/kgrunthal/HD/epta_sim/test_1/")




chainname = 'chain_1'
chain = np.loadtxt("/u/kgrunthal/HD/epta_sim/test_1/" + chainname + '.txt')

burn = int(0.3*chain.shape[0])

#corner.corner(chain[burn:,-4],
#                      bins =30,
#                      plot_datapoints=False, plot_density=True, 
#                      plot_contours=False,fill_contours=False,
#                      show_titles = True, use_math_text=True, verbose=True)

fs = (np.arange(30) + 1) / Tspan
parts = plt.violinplot(chain[burn:,:-4], positions=fs, widths=0.07*fs)
plt.savefig("/u/kgrunthal/HD/epta_sim/test_1/violinplot_single_bin.png")
plt.clf()




ostat = opt_stat.OptimalStatistic(psrs, pta=pta, orf='hd')



chain = np.genfromtxt("/u/kgrunthal/HD/epta_sim/test_1/chain_1.txt")
chain_r = chain[:, :-4]

param_dict = {}
param_dict[pta.param_names[0][:-2]] = chain_r[100]

for i,p in enumerate(pta.param_names):
    param_dict[p] = chain_r[100][i]


print(param_dict)
snr_list = []
freq_list = 10**np.arange(-9, -6, 0.1)
for j in range(len(freq_list)):
    xi, rho, sig, OS, OS_sig = ostat.compute_os(params = param_dict, psd = "spectrum", fgw=freq_list[j])
    print(OS, OS_sig)
    snr_list = np.append(snr_list, OS/OS_sig)

_, _, _, OSpl, OSpl_sig = ostat.compute_os(params = param_dict)
print('powerlaw: A^2 = {} +/- {} with S/N = {}'.format(OSpl, OSpl_sig, OSpl/OSpl_sig))



fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True)
plt.subplots_adjust(hspace=0)

# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
axs[0].set_xscale('log')
axs[0].errorbar(freq_list, snr_list, fmt='ko', ls='')
axs[1].errorbar(1, OSpl/OSpl_sig)
axs[1].set_xticklabels('PL S/N')
plt.savefig("/u/kgrunthal/HD/epta_sim/test_1/snr.png")
plt.clf()


dict(zip([pta.param_names[0][:-2]] * len(pta.param_names), chain_r[100]))

param_dict = {}
param_dict[pta.param_names[0][:-2]] = chain_r[100]


plt.plot(chain_r[100])
plt.savefig("/u/kgrunthal/HD/epta_sim/test_1/chain.png")
plt.clf()
