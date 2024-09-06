# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 16:26:59 2024

@author: kgrun
"""


import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord







Npsr = 100



# galactic longitude
l_rad = rng.uniform(0,2*np.pi, Npsr)

# scale height distribution
zavg = 0.5
zabs = rng.exponential(zavg, Npsr)
z = zabs.copy()
z[1::2] *= -1


plt.hist(z, 10, density=True)
plt.xlabel('scale height')
plt.show()



# heliocentric distance
loc, scale = 1.0, 1.5
dist = np.zeros(Npsr)
for i in range(Npsr):
    dtmp = np.abs(rng.laplace(loc, scale))
    while np.abs(z[i]/dtmp) > 1.:
        dtmp = np.abs(rng.laplace(loc, scale))
    dist[i]=dtmp
        

plt.hist(dist, 100, density=True)
#plt.xlim(0,10)
plt.xlabel('helocentric distance')
plt.show()


# galactic latitude
b_rad = np.arcsin(z/dist)



# convert to RA and DEC


c_gal = SkyCoord(l = l_rad, b=b_rad, unit='rad', frame='galactic')
c = c_gal.icrs
ra, dec = np.zeros(Npsr), np.zeros(Npsr)
for i in range(Npsr):
    if c[i].ra.radian > np.pi:
        ra[i], dec[i] = 2*np.pi-c[i].ra.radian, c[i].dec.radian
    else:
        ra[i], dec[i] = -1*c[i].ra.radian, c[i].dec.radian
        
plt.figure(figsize=(7,9))
plt.subplot(projection="mollweide")
plt.grid(which='both')
plt.errorbar(ra, dec, fmt='ko', ls ='')