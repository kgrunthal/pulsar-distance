#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 13:40:33 2023

@author: kgrunthal
"""

import sys
import pickle, json
import subprocess
import os, glob, json 
import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt

from astropy import units as u

from astropy.coordinates import SkyCoord

import scipy.constants as sc


plt.rcParams['font.family'] = "serif"
plt.rcParams['font.sans-serif'] = "Times"

plt.rcParams['text.usetex']= False
plt.rcParams['xtick.labelsize'] = 11.0
plt.rcParams['ytick.labelsize'] = 11.0
plt.rcParams['axes.labelsize'] = 14.0



SOLAR2S = sc.G / sc.c**3 * 1.98855e30
KPC2S = sc.parsec / sc.c * 1e3
MPC2S = sc.parsec / sc.c * 1e6


def generate_galactic_pulsars(Npsr, datadir, plots=True):
    # galactic longitude
    l_rad = rng.uniform(0,2*np.pi, Npsr)

    # scale height distribution
    zavg = 0.5
    zabs = rng.exponential(zavg, Npsr)
    z = zabs.copy()
    z[1::2] *= -1
    
    # heliocentric distance
    loc, scale = 1.0, 1.5
    dist = np.zeros(Npsr)
    for i in range(Npsr):
        dtmp = rng.laplace(loc, scale)
        while np.abs(z[i]/dtmp) > 1.:
            dtmp = rng.laplace(loc, scale)
        dist[i]=dtmp
    
    # galactic latitude
    b_rad = np.arcsin(z/dist)

    # convert to RA and DEC
    c_gal = SkyCoord(l = l_rad, b=b_rad, unit='rad', frame='galactic')
    c = c_gal.icrs
    ra, dec = c.ra.radian, c.dec.radian
    
    dist_dict = {}
    for ii, ci in enumerate(c):
        cstr = ci.to_string('hmsdms')
        #print cstr
        RAJ = cstr.split(" ")[0].replace("h",":").replace("m",":")[:-1]
        DECJ = cstr.split(" ")[1].replace("d",":").replace("m",":")[:-1]
        cstr = cstr.replace(" ","")
        name = "J"+RAJ[0:2]+RAJ[3:5]+DECJ[0]+DECJ[1:3]+DECJ[4:6]
        
        dist_dict[name] = dist[ii]
    
    with open(datadir + 'distances.json', 'w') as distfile:
        json.dump(dist_dict, distfile, indent=4)
    distfile.close()
    
    if plots == True:
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10,6))
        
        axs[0].hist(z, bins=20)
        axs[1].hist(dist, bins=20)
        
        axs[0].set_xlabel('scale height $z$')
        axs[1].set_xlabel('helocentric distance / kpc')
        plt.savefig(datadir + 'histograms.png', bbox_inches='tight', dpi=400)
        
        plt.clf() 
        ra_plt, dec_plt = np.zeros(Npsr), np.zeros(Npsr)
        for i in range(Npsr):
            if c[i].ra.radian > np.pi:
                ra_plt[i], dec_plt[i] = 2*np.pi-c[i].ra.radian, c[i].dec.radian
            else:
                ra_plt[i], dec_plt[i] = -1*c[i].ra.radian, c[i].dec.radian
            
        plt.figure(figsize=(7,9))
        plt.subplot(projection="mollweide")
        plt.grid(which='both')
        plt.errorbar(ra_plt, dec_plt, fmt='ko', ls ='')
        
        plt.xticks(np.linspace(-np.pi, np.pi, 9), 
                   ['12h', '9h', '6h', '3h', '0h', '21h', '18h', '15h', ''])
        plt.savefig(datadir + 'skydistribution.png', bbox_inches='tight', dpi=400)
        plt.clf()

    return ra, dec





    

def make_fake_pulsar(phi, theta, DIR=""):
    '''
    Makes a fake pulsar par file
    '''
    
    c = SkyCoord(phi,theta,frame='icrs',unit='rad')
    cstr = c.to_string('hmsdms')
    #print cstr
    RAJ = cstr.split(" ")[0].replace("h",":").replace("m",":")[:-1]
    DECJ = cstr.split(" ")[1].replace("d",":").replace("m",":")[:-1]
    cstr = cstr.replace(" ","")
    name = "J"+RAJ[0:2]+RAJ[3:5]+DECJ[0]+DECJ[1:3]+DECJ[4:6]
    print(name)
    output = "PSR           %s\n"%name
    
    output += "RAJ           %s 1\n"%RAJ
    output += "DECJ          %s 1\n"%DECJ
    
    period = 0.001*np.random.uniform(1,5) #seconds
    output += "F0            %0.10f 1\n"%(1.0/period)
    
    output += "PEPOCH        55000.0\n"    
    output += "POSEPOCH      55000.0\n"
    
    dist = np.random.uniform(0.1,3) #kpc
    #output += "PX            %0.5f 1\n"%(1.0/dist)
               
    output += "EPHEM         DE440\n"
    
    output += "CLK           TT(BIPM2021)\n"
    output += "MODE 1\n"
    
    
    filename = "%s%s.par"%(DIR,name)
    with open(filename,'w') as FILE:
        FILE.write(output)

    return filename.encode('ascii','ignore')





def make_parfiles(Npsr, distribution='isotropic', datadir=''):
    # create the par files for the PTA with
    # Npsr - number of pulsars
    # isotropic - if True, distribute pulsars isotropically on the sphere
    if distribution == 'isotropic':
        # Fibonacci sequence on sphere
        i = np.arange(0, Npsr, dtype=float) + 0.5
        golden_ratio = (1 + 5**0.5)/2
        costhetas = 1 - 2*i/Npsr
        thetas = np.arccos(costhetas) - np.pi/2
        phis = np.mod(2 * np.pi * i / golden_ratio, 2*np.pi)
        
        
    elif distribution == 'random':
        costhetas = np.random.uniform(-1., 1., size=Npsr)
        thetas = np.arccos(costhetas) - np.pi/2.
        phis = np.random.uniform(0., 2.*np.pi, size=Npsr)
        
    elif distribution == 'ring':
        phis = np.linspace(0, 2*np.pi, Npsr, endpoint=False)
        thetas = np.ones(Npsr)*(np.pi/2. - np.pi/2.)
        
    elif distribution == 'galactic':
        phis, thetas = generate_galactic_pulsars(Npsr, datadir, plots=True)
        
    for i,n in enumerate(np.arange(Npsr)):
        make_fake_pulsar(phis[i], thetas[i], DIR=datadir+'par/')
    
    return None




##############################################################################
##############################################################################

datadir = "./par/SKA/"



# 1. create pulsars
subprocess.run("mkdir {}".format(datadir+'par').split(' '))

make_parfiles(int(sys.argv[1]), distribution='galactic', datadir=datadir)

