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


from astropy.coordinates import SkyCoord

import scipy.constants as sc


SOLAR2S = sc.G / sc.c**3 * 1.98855e30
KPC2S = sc.parsec / sc.c * 1e3
MPC2S = sc.parsec / sc.c * 1e6


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
        
        
    for i,n in enumerate(np.arange(Npsr)):
        make_fake_pulsar(phis[i], thetas[i], DIR=datadir+'par/')
    
    return None




##############################################################################
##############################################################################

datadir = "./"



# 1. create pulsars
subprocess.run("mkdir {}".format(datadir+'par').split(' '))

make_parfiles(int(sys.argv[1]), distribution='ring', datadir=datadir)

