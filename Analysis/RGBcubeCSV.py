# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 10:57:12 2018

@author: Julia Kent

This code is modeled after Logan's Wright's RGBgeneration.py code for HICO files
but adapted for the HySICS csv img files
"""

#import modules
import csv
import os
import numpy as np
from skimage import exposure
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#Function for creating the datacube and a linear stretched rgb image
def datacube(path,plotting):
    files = [os.path.join(path,f) for f in os.listdir(path) if f.startswith('img')]
    testfiles = files[0:1050]
    
    #Radiance Data
    data_cube = []
    for file in testfiles:
        imgxxxx = open(file)
        reader = csv.reader(imgxxxx)
        data = []
        for row in reader:
            newrow = []
            for i in np.arange(len(row)):
                element = float(row[i])
                newrow.append(element)
            data.append(newrow) 
        imgxxxx.close()
        data_cube.append(data)
        
    data_cube = np.array(data_cube) # (1050,480,640)
    data_cube = np.swapaxes(data_cube,0,1) #(480,1050,640)
       
    #Wavelengh Data
    WLfile = open(os.path.join(path,'wavelength.txt'))
    WLreader = csv.reader(WLfile)
    WL = []
    for row in WLreader:
        newrow = []
        for i in np.arange(len(row)):
            element = float(row[i])
            newrow.append(element)
        WL.append(newrow) 
    WLfile.close()
    WL = np.array(WL)
    
    #Find RGB indices
    Ridx = np.argmin(np.abs(670 - WL))
    Gidx = np.argmin(np.abs(555 - WL))
    Bidx = np.argmin(np.abs(443 - WL))
    
    #Apply a 2% Linear Stretch
    R = data_cube[:,:,Ridx]/np.amax(data_cube[:,:,Ridx])
    p2, p98 = np.percentile(R, (2, 98))  # Makes histogram and finds the value 2% of points in from either end      
    Rscl = exposure.rescale_intensity(R, in_range=(p2, p98)) # Squishes the data in the new smaller range
    
    G = data_cube[:,:,Gidx]/np.amax(data_cube[:,:,Gidx])
    p2, p98 = np.percentile(G,(2, 98))
    Gscl = exposure.rescale_intensity(G, in_range=(p2, p98))
    
    B = data_cube[:,:,Bidx]/np.amax(data_cube[:,:,Bidx])
    p2, p98 = np.percentile(B, (2, 98))
    Bscl = exposure.rescale_intensity(B, in_range=(p2, p98))
    
    #Recombine the stretched rgb channels
    rgb = np.dstack((Rscl,Gscl,Bscl))
    
    #Display RGB image
    if (plotting == True):
        plt.imshow(rgb)
    
    return(data_cube,WL,rgb);

path = '..\LASP py\HySICS Data\Thick_clouds1\img'
#path = '..\LASP py\HySICS Data\Desert_vegetation_clouds'
(data_cube,WL,rgb) = datacube(path, False)

#Modified z-score method for removing outliers 
#This is done because the hyperspectral displays appeared monochrome
def modified_z_score(data):
    threshold = 2.5
    median = np.median(data)
    mad = np.median([np.abs(i - median) for i in data]) #median absolute deviation
    mod_z_scores = [0.6745 * (i - median) / mad for i in data]
    return np.where(np.abs(mod_z_scores) > threshold)

#Creating the Hyperspectral Data Cube
#Makes a hypercube by treating it as a 3d plot with points on 3 surfaces
#More efficient methods may exist
def displaycube(data_cube,rgb):
    imdata = np.flipud(rgb)
    dims = data_cube.shape
    res = 1
    
    #   Define a grid to be plotted on   
    Xi, Zi = np.meshgrid(0.25*np.arange(dims[1]),np.arange(dims[0])) #Front Face (RGB)
    Yi = np.full((dims[0],dims[1]),0,dtype = 'int')
    
    X, Y = np.meshgrid(0.25*np.arange(dims[1]),np.arange(dims[2])) #Top Hyperspectral Face
    Z = np.full((dims[2],dims[1]),dims[0],dtype = 'int')
    
    Z2, Y2 = np.meshgrid(np.arange(dims[0]),np.arange(dims[2])) #Right Hyperspectral Face
    X2 = np.full((dims[2],dims[0]),0.25*dims[1],dtype = 'int')
    
    #Scale the Colors for the Hyperspectral Faces
    C = np.fliplr(np.transpose(data_cube[0,:,:])) #Top Face -- using [-2,:,:] bc [-1,:,:] was all 0s, using[0,:,:] because rgb image has been flipud
    outliers_top = modified_z_score(C)
    C[outliers_top] = False #Mask out extreme values
    
    C2 = np.transpose(data_cube[:,-1,:])  #Right Face
    outliers_right = modified_z_score(C2)
    C2[outliers_right] = False
    C2 = np.fliplr(C2)
    
    CCmin = np.max([np.min(C),np.min(C2)]) #Bounds for Normalizing Color Map
    CCmax = np.min([np.max(C),np.max(C2)])
    Cnorm = mpl.colors.Normalize(vmin = CCmin, vmax = CCmax) #Normalize color map
    
    #Finally Display the Hyperspectral Data Cube
    hypercube = plt.figure(figsize = (2,5))
    ax = hypercube.gca(projection = '3d',aspect = 'equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    #Fill the surfaces
    ax.plot_surface(Xi, Yi, Zi, rstride = res, cstride = res, facecolors = imdata, linewidth=0, shade = False, antialiased = False) #Front Face    
    ax.plot_surface(X, Y, Z, rstride = res, cstride = res, facecolors = mpl.cm.jet(Cnorm(C)), linewidth=0) #Top
    ax.plot_surface(X2, Y2, Z2, rstride = res, cstride = res, facecolors = mpl.cm.jet(Cnorm(C2)), linewidth=0) #Right
    
    ax.azim = 285
    ax.elev = 15
    scaling = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz']);
    ax.auto_scale_xyz(*[[np.min(scaling), np.max(scaling)]]*3)
    plt.axis('off')
    return;
    
displaycube(data_cube, rgb)