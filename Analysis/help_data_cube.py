import csv
import scipy
import numpy as np
from skimage import exposure
import matplotlib as mpl
import matplotlib.pyplot as plt

#Radiance Data
def find_data_cube_radiances(testfiles):
    data_cube = []
    for file in testfiles:
        imgx = open(file)
        reader = csv.reader(imgx)
        data = []
        for row in reader:
            newrow = []
            for i in np.arange(len(row)):
                element = float(row[i])
                newrow.append(element)
            data.append(newrow) 
        imgx.close()
        data_cube.append(data)
        
    data_cube = np.array(data_cube) # (1050,480,640)
    data_cube = np.swapaxes(data_cube,0,1) #(480,1050,640)
    return data_cube;


#Wavelengh Data
def read_HySICS_wl(phase_dict):
    hysics_wl = scipy.io.readsav(phase_dict['HySICS_wl_path'])
    hysics_wl = hysics_wl['wlsample']
    hysics_wl = hysics_wl[0:-27] if phase_dict['phase']=='Liquid Water' else hysics_wl[35:-27]
    hysics_wl = [x+4 for x in hysics_wl]
    return hysics_wl;


def two_percent_linear_stretch(x,idx):
    c = x[:,:,idx]/np.amax(x[:,:,idx])
    p2, p98 = np.percentile(c, (2, 98))
    c_scaled = exposure.rescale_intensity(c, in_range=(p2, p98))
    return c_scaled;


def rgb_idx(wl_list):
    r_idx = np.argmin([np.abs(670 - x) for x in wl_list])
    g_idx = np.argmin([np.abs(555 - x) for x in wl_list])
    b_idx = np.argmin([np.abs(443 - x) for x in wl_list])
    rgb_idx_list = [r_idx,g_idx,b_idx]
    return (rgb_idx_list); 


def apply_stretch(data_cube,wl_list):
    rgb_idx_list = rgb_idx(wl_list)
    
    red = two_percent_linear_stretch(data_cube,rgb_idx_list[0])
    green = two_percent_linear_stretch(data_cube,rgb_idx_list[1])
    blue = two_percent_linear_stretch(data_cube,rgb_idx_list[2])

    rgb = np.dstack((red,green,blue))
    return rgb;    


def modified_z_score(data): #This is done because the hyperspectral displays appeared monochrome
    threshold = 2.5
    median = np.median(data)
    mad = np.median([np.abs(i - median) for i in data]) #median absolute deviation
    mod_z_scores = [0.6745 * (i - median) / mad for i in data]
    return np.where(np.abs(mod_z_scores) > threshold)


def displaycube(data_cube,rgb): #Makes a hypercube by treating it as a 3d plot with points on 3 surfaces
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
    
    #Display the Hyperspectral Data Cube
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