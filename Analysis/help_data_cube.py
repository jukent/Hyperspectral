import csv
import numpy as np
from skimage import exposure
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
def find_wl_data(path):
    wl_file = open(os.path.join(path,'wavelength.txt'))
    wl_reader = csv.reader(wl_file)
    wl_list = []
    for row in wl_reader:
        newrow = []
        for i in np.arange(len(row)):
            element = float(row[i])
            newrow.append(element)
        wl_list.append(newrow) 
    WLfile.close()
    wl_list = np.array(wl_list)
    
    #HySICS 4nm wavelength shift
    wl_list = [x+4 for x in wl_list]
    #Filter wl_list but by how much????

    #Find RGB indices
    r_idx = np.argmin(np.abs(670 - wl_list))
    g_idx = np.argmin(np.abs(555 - wl_list))
    b_idx = np.argmin(np.abs(443 - wl_list))
    rgb_idx_dict = {'red_index':r_idx,'green_index':g_idx,'blue_index':b_idx}
    return (wl_list,rgb_idx_dict);

def two_percent_linear_stretch(x,idx):
    c = x[:,:,idx]/np.amax(x[:,:,idx])
    p2, p98 = np.percentile(c, (2, 98))
    c_scaled = exposure.rescale_intensity(c, in_range=(p2, p98))
    return c_scaled;

#Apply a 2% Linear Stretch
def apply_stretch(data_cube,rgb_idx)
    red = two_percent_linear_stretch(data_cube,rgb_idx['red_index'])
    green = two_percent_linear_stretch(data_cube,rgb_idx['green_index'])
    blue = two_percent_linear_stretch(data_cube,rgb_idx['blue_index'])

    #Recombine the stretched rgb channels
    rgb = np.dstack((red,green,blue))
    return rgb;    


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