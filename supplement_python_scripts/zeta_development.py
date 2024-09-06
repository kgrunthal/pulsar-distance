# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 16:26:58 2024

@author: kgrun
"""


import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt


def zeta(cosMu, pd):
    KPC2S = sc.parsec / sc.c * 1e3
    SOLAR2S = sc.G / sc.c**3 * 1.98855e30

    w0 = 2*np.pi*1e-8
    mc = 10**(8.5)
    mc *= SOLAR2S
    pd *= KPC2S
    
    fac1 = 256 / 5 * mc ** (5 / 3) * w0 ** (8 / 3)
    return (1 + fac1 * pd * (1 - cosMu)) ** (-3 / 8)

cosMu = np.linspace(-1., 1, 10)
pd_1 = 1.*np.ones(10)
pd_2 = 2.*np.ones(10)
pd_3 = 3.*np.ones(10)


plt.errorbar(cosMu, zeta(cosMu, pd_1), ls='', fmt='kx')
plt.errorbar(cosMu, zeta(cosMu, pd_2), ls='', fmt='rx')
plt.errorbar(cosMu, zeta(cosMu, pd_3), ls='', fmt='bx')
plt.show()


