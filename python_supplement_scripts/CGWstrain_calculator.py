# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 15:42:23 2024

@author: kgrun
"""


import numpy as np
import scipy.constants as sc

SOLAR2S = sc.G / sc.c**3 *1.98855e30
MPC2S = sc.parsec / sc.c * 1e6

def strain(fgw, logMc, dlum):
    Mc = (10**logMc)*SOLAR2S
    dlum = dlum*MPC2S
    h = 2 * Mc**(5/3) * (np.pi*fgw)**(2/3) /dlum
    return h

def scaletodlum(fgw, logMc, dlum=15):
    Mc = (10**logMc)*SOLAR2S
    dlum = dlum*MPC2S
    h = 2 * Mc**(5/3) * (np.pi*fgw)**(2/3) /dlum
    return h


def kappa(lmc, lmc_ref = 8.7):
    logkappa = (5/3)*(lmc - lmc_ref)
    kappa = 10**logkappa
    print('kappa check:', lmc, lmc_ref, logkappa, kappa)
    return kappa

print('strain 8.7, d = 15 Mpc ', strain(22.3e-9, 8.7, 15))

print()
print(5*strain(22.3e-9, 8.7, 15), strain(22.3e-9, 9.12, 15))

dl_new = kappa(8.7, 9.119) * 15
print(dl_new)
print(strain(22.3e-9, 8.7, dl_new))

dl_90 = kappa(9.0, 8.7) * 15

print(dl_90)
print(strain(22.3e-9, 9.0, dl_90))


print('zhu ', strain(1e-8, np.log10(8.77*1e8), 16.5))