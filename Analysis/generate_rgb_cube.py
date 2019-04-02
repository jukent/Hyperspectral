# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 10:57:12 2018

@author: Julia Kent

This code is modeled after Logan's Wright's RGBgeneration.py code for HICO files
but adapted for the HySICS csv img files
"""

#import modules
import os
import help_data_cube as dc


#Function for creating the datacube and a linear stretched rgb image
def generate_datacube(path):
    files = [os.path.join(path,f) for f in os.listdir(path) if f.startswith('img')]
    testfiles = files[0:1050]
    
    data_cube = dc.find_data_cube_radiances(testfiles)
    wl_list, rgb_idx = dc.find_wl_data(path)  
    rgb = dc.apply_stretch(data_cube,rgb_idx)
    
    return(data_cube,wl_list,rgb);

path = '..\LASP py\HySICS Data\Thick_clouds1\img'
#path = '..\LASP py\HySICS Data\Desert_vegetation_clouds'
(data_cube,wl_list,rgb_idx,rgb) = datacube(path)

plt.imshow(rgb)    
dc.displaycube(data_cube, rgb)
