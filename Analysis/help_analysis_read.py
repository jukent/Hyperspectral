# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 02:27:06 2019

@author: julia
"""

import numpy as np
import csv
import scipy.io
import os

def read_HySICS(phase_dict):
    wlHy = scipy.io.readsav(phase_dict['HySICS_wl_path'])
    wlHy = wlHy['wlsample']
    wlHyf=wlHy[0:-27] if phase_dict['phase']=='Liquid Water' else wlHy[35:-27]
    cloud = np.load(phase_dict['HySICS_data_path'])
    cloudf=cloud[0:-27] if phase_dict['phase'] =='Liquid Water' else cloud[35:-27]
    hysics_dict = {'wl':wlHyf,'data':cloudf,'phase':phase_dict['phase']}
    return hysics_dict;


def read_HySICS_wl(phase_dict):
    hysics_wl = scipy.io.readsav(phase_dict['HySICS_wl_path'])
    hysics_wl = hysics_wl['wlsample']
    hysics_wl = hysics_wl[0:-27] if phase_dict['phase']=='Liquid Water' else hysics_wl[35:-27]
    return hysics_wl;


def read_LRT(file,phase):
    data = [x for x in csv.reader(open(file,'r'),delimiter='\t')]  
    data = data[5:]
    wl=[]           #Wavelength (nm)
    dirL=[]         #Direct Radiance (W/m^2/nm/sr)
    
    loops = int(len(data)/3)
    for n in range(loops):
        if phase =='Liquid Water':           
            nwl = np.float(data[3*n+1][0][1:9])
            ndirL = np.float(data[3*n+3][0][9:24]) 
        else:
            nwl = np.float(data[3*n][0][1:9])
            ndirL = np.float(data[3*n+2][0][9:24])

        wl.append(nwl)
        dirL.append(ndirL)

    dirLW = [i*10**-3 for i in dirL] #convert to W/m^2/nm/sr
    return(wl,dirLW);
    

def read_all_LRT(phase_dict):
    LRTfiles = os.listdir(phase_dict['LRT_path']) 
    LRT_vals = []
    for f in LRTfiles:
        if f.endswith(".dat"):
            (wl,LRT_val) = read_LRT(phase_dict['LRT_path']+f,phase_dict['phase'])
            LRT_vals.append(LRT_val)  
    LRT_dict = {'wl':wl,'data':LRT_vals,'phase':phase_dict['phase']}
    return LRT_dict;

    
def read_verbose(file,phase):
    data = [x for x in open(file,'r')]
    a = data.index('Using new intensity correction, with phase functions\n')
    n=4 if phase=='Liquid Water' else 5
    tau = data[a:][67].split('|')[n].split()[0]
    return(tau);


def read_solar(file):
    d = np.loadtxt(file, delimiter="\t")
    solar_wl=d[:,0]
    solar_flux=d[:,1]
    solar_flux = np.array([i/1000 for i in solar_flux]) 
    solar_dict = {'wl':solar_wl,'flux':solar_flux}
    return solar_dict;