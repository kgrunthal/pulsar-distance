# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 08:53:10 2024

@author: kgrun
"""


import numpy as np
import matplotlib.pyplot as plt

'''
data = np.loadtxt('GWBspectrum_singlebin.txt')

plt.plot(np.log10(data[0]),data[1], ls='', marker='o')
#plt.xlim(-8.0, -7.9)
#plt.ylim(-15.2, -14.8)
plt.show()
'''

howml = 2.
dur = 12*365.25*24*3600


def broken(f, Amp, gamma, beta, f0, power):
    f1yr = 1 / 3.16e7
    alpha = -0.5 * (gamma - 3)
    hcf = Amp * (f / f1yr) ** (alpha)
    si = alpha - beta
    print(si)

    hcf /= (1 + (f / f0) ** (power * si)) ** (1 / power)
    C = 1 / 96 / np.pi**2 * hcf**2 / f**3 * dur * howml
    return hcf


def powerlaw(f, Amp, gamma):
    f1yr = 1 / 3.16e7
    alpha = -0.5 * (gamma - 3)
    
    hcf = Amp * (f / f1yr) ** (alpha)
    C = 1 / 96 / np.pi**2 * hcf**2 / f**3 * dur * howml
    return hcf


def singlebin(f, b):
    spec = 1e-50*np.ones(len(f))
    spec[b] = 6e-15*np.ones(1)
    
    return spec


fmin = 1/(10*365.25*24*3600)

ftest = np.logspace(np.log10(fmin), np.log10(20*fmin), 20)

Amp_normal = 1e-14
Amp_broken = 1e-14
gamma = 13/3
beta = 13/3
f0 = 1.2e-8

plt.plot(np.log10(ftest), broken(ftest, Amp_broken, gamma, beta, f0, 2), ls='', marker='x')
plt.plot(np.log10(ftest), powerlaw(ftest, Amp_normal, gamma), ls='', marker='o')
plt.plot(np.log10(ftest), singlebin(ftest, 19), ls='', marker='o')
plt.ylim(1e-15, 1e-13)

plt.xlabel('log$_{10}f$', fontsize=14)
plt.ylabel('$h_c(f)$',  fontsize=14)
plt.yscale('log')
plt.show()
    
    

