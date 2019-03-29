# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 02:27:06 2019

@author: julia
"""

import numpy as np
import csv
import scipy.io
import os

def read_HySICS(file):
    wlHy = scipy.io.readsav(file)
    wlHy = wlHy['wlsample']
    wlHyf=wlHy[0:-27]
    thincloud = np.load('thincloud.npy')
    thincloudf=thincloud[0:-27]
    hysics_dict = {'wl':wlHyf,'data':thincloudf}
    return hysics_dict;


def read_LRT(file):
    data = [x for x in csv.reader(open(file,'r'),delimiter='\t')]  
    data = data[5:]
    wl=[]           #Wavelength (nm)
    dirL=[]         #Direct Radiance (W/m^2/nm/sr)
    
    loops = int(len(data)/3)
    for n in range(loops):
        nwl = np.float(data[3*n+1][0][1:9])
        ndirL = np.float(data[3*n+3][0][9:24])      
        wl.append(nwl)
        dirL.append(ndirL)    
    dirLW = [i*10**-3 for i in dirL] #convert to W/m^2/nm/sr
    return(wl,dirLW);
    

def read_all_LRT(LRT_path):
    LRTfiles = os.listdir(LRT_path) 
    LRT_vals = []
    for f in LRTfiles:
        if f.endswith(".dat"):
            (wl,LRT_val) = read_LRT(LRT_path+'/'+f)
            LRT_vals.append(LRT_val)  
    LRT_dict = {'wl':wl,'data':LRT_vals}
    return LRT_dict;

    
def read_verbose(file):
    data = [x for x in open(file,'r')]
    a = data.index('Using new intensity correction, with phase functions\n')
    tau = data[a:][67].split('|')[4].split()[0]
    return(tau);


def read_solar(file):
    d = np.loadtxt(file, delimiter="\t")
    solar_wl=d[:,0]
    solar_flux=d[:,1]
    solar_flux = np.array([i/1000 for i in solar_flux]) 
    solar_dict = {'wl':solar_wl,'flux':solar_flux}
    return solar_dict;