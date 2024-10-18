#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 15:20:28 2023

@author: kgrunthal
"""

import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

SOLAR2S = sc.G / sc.c**3 * 1.98855e30
KPC2S = sc.parsec / sc.c * 1e3
MPC2S = sc.parsec / sc.c * 1e6



def omega_p_noevo(pd, mc=1e10, fgw=1e-8 , cosMu=0.):
    mc *= SOLAR2S  # convert from solar masses to seconds
    pd *= KPC2S

    w0 = np.pi * fgw
    fac1 = 256 / 5 * mc ** (5 / 3) * w0 ** (8 / 3)
    omega_p = w0 * (1 + fac1 * pd * (1 - cosMu)) ** (-3 / 8)
    
    return omega_p-w0


def omega_p_evo(toas, pd,  mc=10**9., fgw=5e-8 , cosMu=0.):
    mc *= SOLAR2S  # convert from solar masses to seconds
    pd *= KPC2S
    toas = toas * 86400 
    tp = toas - pd * (1 - cosMu)
    w0 = np.pi * fgw
    fac1 = 256 / 5 * mc ** (5 / 3) * w0 ** (8 / 3)
    
    omega = w0 * (1 - fac1 * toas) ** (-3 / 8)
    omega_p = w0 * (1 - fac1 * tp) ** (-3 / 8)
    
    return omega-omega_p

distances = np.arange(0.5,2.5,0.5)

toas = np.arange(50000,53652,14)

for d in distances:
    plt.plot(toas, omega_p_evo(toas, d), label='{} kpc'.format(d))
    
plt.xlabel('TOA')
plt.ylabel('$\omega_p - \omega_0$')
plt.title('Mc = 1e8, fgw = {}'.format(5e-8))
plt.legend(loc='best')
plt.show()


