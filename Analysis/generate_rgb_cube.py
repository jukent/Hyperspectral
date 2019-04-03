# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 10:57:12 2018

@author: Julia Kent

This code is modeled after Logan's Wright's RGBgeneration.py code for HICO files
but adapted for the HySICS csv img files
"""

#import modules
import os
import matplotlib.pyplot as plt
import help_data_cube as dc

HySICS_wl_path = '../../../CU Boulder/LASP/LASP py/HySICS Data/WLHysics.sav'
HySICS_LW_data_path = '..\..\LASP py\HySICS Data\Thick_clouds1\img'
HySICS_IC_data_path = '..\LASP py\HySICS Data\Desert_vegetation_clouds'

water_dict = {'HySICS_wl_path':HySICS_wl_path,'HySICS_data_path':HySICS_LW_data_path,'phase':'Liquid Water'}
ice_dict = {'HySICS_wl_path':HySICS_wl_path,'HySICS_data_path':HySICS_IC_data_path,'phase':'Ice'}

#phase_dict=water_dict
phase_dict=ice_dict

def generate_datacube(phase_dict):
    path = phase_dict['HySICS_data_path']
    files = [os.path.join(path,f) for f in os.listdir(path) if f.startswith('img')]
    testfiles = files[0:1050]
    
    data_cube = dc.find_data_cube_radiances(testfiles)
    wl_list = dc.read_HySICS_wl(phase_dict)
    rgb = dc.apply_stretch(data_cube,wl_list)
    
    return(data_cube,wl_list,rgb);


(data_cube,wl_list,rgb) = generate_datacube(phase_dict)

plt.imshow(rgb)    
dc.displaycube(data_cube, rgb)
