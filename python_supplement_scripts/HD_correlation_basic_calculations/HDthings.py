#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 14:49:05 2023

@author: kgrunthal
"""

import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt
import astropy.coordinates as coord

SOLAR2S = sc.G / sc.c**3 * 1.98855e30
KPC2S = sc.parsec / sc.c * 1e3
MPC2S = sc.parsec / sc.c * 1e6



def HD(theta):
    x = (1-np.cos(theta))/2 
    t1 = x * np.log(x)
    t2 = x/6.
    return (3/2.) * (t1 - t2 +1/3.)

def angular_dist(phi1, phi2, theta1, theta2):
    
    ra1 = phi1
    ra2 = phi2
    dec1 = np.pi/2 - theta1
    dec2 = np.pi/2 - theta2
    
    theta = np.arccos(np.sin(dec1)*np.sin(dec2)+np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2))
   
    return theta
    


def CGW(toas,
        ptheta, pphi,
        gwtheta, gwphi,
        mc, dist, fgw, phase0,
        psi, inc,
        zeta=None, #relative difference between omega and omega_p
        pdist=1.0,
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
    #print("omhat = ", omhat)

    # various factors invloving GW parameters
    fac1 = 256 / 5 * mc ** (5 / 3) * w0 ** (8 / 3)
    fac2 = 1 / 32 / mc ** (5 / 3)
    fac3 = mc ** (5 / 3) / dist

    # pulsar location

    # use definition from Sesana et al 2010 and Ellis et al 2012
    phat = np.array([np.sin(ptheta) * np.cos(pphi), np.sin(ptheta) * np.sin(pphi), np.cos(ptheta)])
    #print("phat = ", phat)
    fplus = 0.5 * (np.dot(m, phat) ** 2 - np.dot(n, phat) ** 2) / (1 + np.dot(omhat, phat))
    fcross = (np.dot(m, phat) * np.dot(n, phat)) / (1 + np.dot(omhat, phat))
    cosMu = -np.dot(omhat, phat)
    #print(" cosMu = ", cosMu)
    


    # convert units
    pd = pdist * KPC2S  # convert from kpc to seconds
    toas = toas * 86400 - tref
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
            
        else:
            omega_p = np.random.uniform((zeta-0.05)*omega, (zeta+0.05)*omega)


        # phases
        phase = phase0 + omega * toas
        phase_p = np.random.uniform(0,2*np.pi) + omega_p * toas
        #phase_p = phase0 + fac2 * (w053 - omega_p ** (-5 / 3)) + omega_p * toas

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
    alpha = fac3 / omega ** (1 / 3.)
    alpha_p = fac3 / omega_p ** (1 / 3.)

    # define rplus and rcross
    rplus = alpha * (At * cos2psi + Bt * sin2psi)
    rcross = alpha * (-At * sin2psi + Bt * cos2psi)
    rplus_p = alpha_p * (At_p * cos2psi + Bt_p * sin2psi)
    rcross_p = alpha_p * (-At_p * sin2psi + Bt_p * cos2psi)

    # residuals
    if psrTerm:
        res = fplus * (rplus_p - rplus) + fcross * (rcross_p - rcross)
    else:
        res = -fplus * rplus - fcross * rcross
    
    return res/np.max(res)




def sort_and_bin(xi, rho, nbins):
    # sort the cross-correlations by xi
    idx = np.argsort(xi)

    xi_sorted = xi[idx]
    rho_sorted = rho[idx]
    
    # bin the cross-correlations so that there are the same number of pairs per bin
    npairs = int(len(idx)/nbins)

    xi_mean = []
    xi_err = []

    rho_avg = []

    i = 0
    while i < len(xi_sorted):
    
        xi_mean.append(np.mean(xi_sorted[i:npairs+i]))
        xi_err.append(np.std(xi_sorted[i:npairs+i]))

        r = np.mean(rho_sorted[i:npairs+i])
        rho_avg.append(r)
   
        i += npairs
    
    xi_mean = np.array(xi_mean)
    xi_err = np.array(xi_err)
    
    return xi_mean, xi_err, rho_avg



#########

'''
### One realisation, pulsars in a ring ###

Npulsars = 20

CGWs = []
PHIs = np.linspace(0, 2*np.pi, Npulsars, endpoint=False)

#test
#PHIs = np.random.uniform(0, 2*np.pi, Npulsars)
#COSTHETAs = np.random.uniform(-1,1,Npulsars)
#THETAs = np.arccos(COSTHETAs)

toas = np.arange(50000,53652,14)
distances = np.random.uniform(2.8,3.2, Npulsars)

fig = plt.figure()
ax = fig.add_subplot()
plt.plot(np.sin(np.linspace(0, 2*np.pi, 200, endpoint=False)), np.cos(np.linspace(0, 2*np.pi, 200, endpoint=False)) , ls ='--', color='lightgray')
plt.errorbar(np.sin(PHIs), np.cos(PHIs), fmt='ro')
ax.set_aspect('equal')
plt.show()




for i, pphi in enumerate(PHIs):
    cgw = CGW(toas,
               np.pi/2, pphi,  # ptheta pphi
               np.pi, 0., #CGWtheta CGWphi
               10**8., 15., 1e-7,
               0., 0., 0.,
               #pdist=1.,
               pdist=distances[i],
               evolve=False,
               phase_approx=True,
               psrTerm=True)
    plt.plot(toas, cgw, ls='-', marker='x', label=pphi)
    
    CGWs.append(cgw)

plt.show()



separations = np.linspace(0, np.pi, int(Npulsars/2)+1)
sepvals = []
corrvals = []
correlations = np.zeros(int(Npulsars/2)+1)
contributions = np.zeros(int(Npulsars/2)+1)

binsize = (separations[2]-separations[1])/2.

for i ,wave in enumerate(CGWs):
    n = i
    while n < len(CGWs):
        
        #asep = min( PHIs[n]-PHIs[i], 2*np.pi-(PHIs[n]-PHIs[i]) )
        
        asep = angular_dist(PHIs[i], PHIs[n], np.pi/2., np.pi/2)
 
        
        sepvals.append(asep)
        corrvals.append(np.correlate(wave, CGWs[n]))
        
        mask = np.argwhere((separations > (asep-binsize)) & (separations < (asep + binsize)) )
        correlations[mask] += np.correlate(wave, CGWs[n])
        contributions[mask] += 1

        n +=1


correlations = correlations/contributions


#plt.errorbar(np.array(sepvals)*180/np.pi, corrvals, ls='', fmt='kx')

plt.errorbar(separations*180./np.pi, correlations/np.max(correlations), ls ="", fmt='ko')
#plt.plot(np.arange(0,180,1), HD(np.linspace(0.001, np.pi, 180))*2)

plt.title("Number of pulsars = {}".format(Npulsars))
plt.xlabel('angular separation / degrees')
plt.show()
'''





def pulsar_ring_realisations(Npulsars, Nreal, Nbins, pdist, zeta = None, evolve=False, phase_approx=True, psrTerm=True):
    output_realisations = np.zeros((Nreal,Nbins))

    for r in range(Nreal):
        PHIs = np.random.uniform(0, 2*np.pi, Npulsars)
        PHIs = np.linspace(0, 2*np.pi, Npulsars, endpoint=False)
        
        #COSTHETAs = np.random.uniform(-1,1,Npulsars)
        #THETAs = np.arccos(COSTHETAs)
        THETAs = np.ones(Npulsars)*np.pi/2
        
        if len(pdist)==2:
            distances = np.random.uniform(pdist[0],pdist[1], Npulsars)
            
        else:
            distances = np.ones(Npulsars)*pdist[0]
    
        toas = np.arange(50000,53652,7)
        CGWs = []
        
        for i, pphi in enumerate(PHIs):
            cgw = CGW(toas,
                      THETAs[i], pphi,  # ptheta pphi
                      np.pi, 0., #CGWtheta CGWphi
                      10**9., 15., 1e-8,
                      0., 0., 0.,
                      zeta=zeta,
                      pdist=distances[i],
                      evolve=evolve,
                      phase_approx=phase_approx,
                      psrTerm=psrTerm)
        
            CGWs.append(cgw)

              
        indseps = []
        corrs = []
        for i,wave in enumerate(CGWs):
            n = i
            while n < len(CGWs):
                
                asep = angular_dist(PHIs[i], PHIs[n], THETAs[i], THETAs[n])
                indseps.append(asep)
                corrs.append(np.correlate(wave, CGWs[n]))
                n+=1
            
            
        xi_mean, xi_err, rho_avg = sort_and_bin(np.array(indseps), np.array(corrs), Nbins)

        #plt.errorbar(xi_mean*180./np.pi, rho_avg/max(rho_avg)/2, xerr=xi_err, ls='', fmt='kx')
        #plt.show()

  
        output_realisations[r]= rho_avg/max(rho_avg)
        #plt.errorbar(separations*180./np.pi, correlations/np.max(correlations), ls='', marker='o', markersize=5)
    

    #plt.title('{} pulsars, {} realisations, psrterm'.format(Npulsars, Nreal))
    #plt.show()

    mean = np.mean(output_realisations, axis=0)
    std = np.std(output_realisations, axis=0)
    
    return mean, std, xi_mean, xi_err



Npulsars = 20
Nbins= 10
Nreal = 500


mean1, std1, xi_mean, xi_err = pulsar_ring_realisations(Npulsars, Nreal, Nbins, [0.8,1.2], psrTerm=True)
mean2, std2, _, _ = pulsar_ring_realisations(Npulsars, Nreal, Nbins, [4.,6.])
mean3, std3, _, _ = pulsar_ring_realisations(Npulsars, Nreal, Nbins, [8.,12])


plt.errorbar(xi_mean*180./np.pi, mean1, xerr=xi_err, yerr=std1, capsize=3, linewidth=3, ls='', fmt = 'r', label='(1.0 +/- 0.2) kpc')
plt.errorbar(xi_mean*180./np.pi+2., mean2, xerr=xi_err, yerr=std2, capsize=3, linewidth=3, ls='', fmt = 'b', label='(5.0 +/- 1.0) kpc')
plt.errorbar(xi_mean*180./np.pi+4., mean3, xerr=xi_err,yerr=std3, capsize=3, linewidth=3, ls='', fmt = 'g', label='(10.0 +/- 2.0) kpc')

#plt.plot(np.linspace(0.01,np.pi,200)*180./np.pi, HD(np.linspace(0,np.pi,200)))

plt.legend(loc='best')
plt.xlabel('angular separation / degrees')
plt.title('different distances')
#plt.savefig('/fpra/mkat/01/users/kgrunthal/Documents/HD/distances_nopsrterm.png', bbox_inches='tight', dpi=400)
plt.show()



'''
mean1, std1, separations = pulsar_ring_realisations(Npulsars, Nreal, [1.], zeta=0.9)
mean2, std2, _ = pulsar_ring_realisations(Npulsars, Nreal, [1.], zeta=0.6)
mean3, std3, _ = pulsar_ring_realisations(Npulsars, Nreal, [1.], zeta=0.3)


plt.errorbar(separations*180./np.pi, mean1, yerr=std1, capsize=3, linewidth=3, ls='', fmt = 'r', label='$\zeta = 0.9$')
plt.errorbar(separations*180./np.pi+2., mean2, yerr=std2, capsize=3, linewidth=3, ls='', fmt = 'b', label='$\zeta = 0.6$')
plt.errorbar(separations*180./np.pi+4., mean3, yerr=std3, capsize=3, linewidth=3, ls='', fmt = 'g', label='$\zeta = 0.3$')

plt.legend(loc='best')
plt.xlabel('angular separation / degrees')
plt.title('$\omega_p = \zeta \omega_0$')
#plt.savefig('/fpra/mkat/01/users/kgrunthal/Documents/HD/scaled_omegap.png', bbox_inches='tight', dpi=400)
plt.show()




mean1, std1, separations = pulsar_ring_realisations(Npulsars, Nreal, [0.8,1.2], phase_approx=False)
mean2, std2, _ = pulsar_ring_realisations(Npulsars, Nreal, [4.,6.], phase_approx=False)
mean3, std3, _ = pulsar_ring_realisations(Npulsars, Nreal, [8.,12], phase_approx=False)


plt.errorbar(separations*180./np.pi, mean1, yerr=std1, capsize=3, linewidth=3, ls='', fmt = 'r', label='(1.0 +/- 0.2) kpc')
plt.errorbar(separations*180./np.pi+2., mean2, yerr=std2, capsize=3, linewidth=3, ls='', fmt = 'b', label='(5.0 +/- 1.0) kpc')
plt.errorbar(separations*180./np.pi+4., mean3, yerr=std3, capsize=3, linewidth=3, ls='', fmt = 'g', label='(10.0 +/- 2.0) kpc')
plt.legend(loc='best')
plt.xlabel('angular separation / degrees')
plt.title('$\omega_p = \omega_0$')
plt.savefig('/fpra/mkat/01/users/kgrunthal/Documents/HD/distances.png', bbox_inches='tight', dpi=400)
#plt.savefig('/fpra/mkat/01/users/kgrunthal/Documents/HD/sameomega', bbox_inches='tight', dpi=400)
plt.show()

'''
