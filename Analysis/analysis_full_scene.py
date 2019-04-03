import numpy as np
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
from help_analysis_read import read_solar
from help_analysis import calc_HySICS_reflectance, disagreement_all_LRT, mask_ground
import help_data_cube as dc

solar_path = '../../../CU Boulder/LASP/LASP py/Solar Data/Aug14SolarRef.dat'
HySICS_wl_path = '../../../CU Boulder/LASP/LASP py/HySICS Data/WLHysics.sav'
HySICS_LW_data_path = 'D:/Flight2_Corrected_Data_Cubes/Desert_vegetation_clouds/'
HySICS_IC_data_path = 'D:/Flight2_Corrected_Data_Cubes/Thick_clouds1/img/'
LRT_LW_path = '../../../CU Boulder/LASP/LASP py/WC files/Outputs/'
verbose_LW_path = '../../../CU Boulder/LASP/LASP py/WC files/Verbose/'
LRT_IC_path = '../../../CU Boulder/LASP/LASP py/ICfiles/Repeat for Baum/Outputs/ghm/'#solid-column/'#rough-aggregate/'
verbose_IC_path = '../../../CU Boulder/LASP/LASP py/ICfiles/Verbose/ghm/'#solid-column/'#rough-aggregate/'

water_dict = {'solar_path':solar_path,'HySICS_wl_path':HySICS_wl_path,'HySICS_data_path':HySICS_LW_data_path, \
              'LRT_path':LRT_LW_path,'LRT_verbose_path':verbose_LW_path, 'phase':'Liquid Water'}
ice_dict = {'solar_path':solar_path,'HySICS_wl_path':HySICS_wl_path,'HySICS_data_path':HySICS_IC_data_path, \
              'LRT_path':LRT_IC_path,'LRT_verbose_path':verbose_IC_path, 'phase':'Ice'}

#-------------------------
wl_num = 15
phase_dict = water_dict
#phase_dict = ice_dict
#-------------------------

solar_dict = read_solar(solar_path)
hysics_wl = dc.read_HySICS_wl(phase_dict) #values already shifted by 4 nm

path = phase_dict['HySICS_data_path']
files = [os.path.join(path,f) for f in os.listdir(path) if f.startswith('img')]
testfiles = files[0:1050]
data_cube = dc.find_data_cube_radiances(testfiles)
data_cube = np.load('data_cube_water.npy')
data_cube = np.load('data_cube_ice.npy')

xx = len(data_cube[:,0,0])
yy = len(data_cube[0,:,0])

tau_scene = np.zeros((xx,yy), dtype = float)
reff_scene = np.zeros((xx,yy), dtype = float)
sdev_scene = np.zeros((xx,yy), dtype = float)
mask = np.zeros((xx,yy), dtype = float)

for y in tqdm(np.arange(0,yy)):
    for x in np.arange(0,xx):
        pixel = data_cube[x,y,:]
        hysics_data = pixel[0:-27] if phase_dict['phase'] == 'Liquid Water' else pixel[35:-27]
        hysics_dict = {'wl':hysics_wl,'data':pixel,'phase':phase_dict['phase']}
        
        st_dev, mask_pixel  = mask_ground(hysics_dict)
        print(mask_pixel)
        sdev_scene[x,y] = st_dev
        if (mask_pixel == True): #Ground -- No Clouds
            tau_scene[x,y] = 'NaN'
            reff_scene[x,y] = 'NaN'                        
        elif (mask_pixel == False): #No Ground -- Yes Clouds
            hysics_dict= calc_HySICS_reflectance(hysics_dict,solar_dict)
            analysis_dict = disagreement_all_LRT(hysics_dict, wl_num, phase_dict) 
            tau_scene[x,y] = analysis_dict['best_COT']
            reff_scene[x,y] = analysis_dict['best_reff']


##Check Mask
#fig, axs = plt.subplots(1,2,sharex=True)
#ax1 = axs[0]
#ax1.imshow(rgb)
#
#rgb_masked = rgb*mask_ground
#ax2 = axs[1]
#ax2.imshow(rgb_masked)
#
#
##Plot tau_scene and reff_scene -- Needs Testing
#fig, axs = plt.subplots(1,2,sharex=True)
#
#ax1  = axs[0]
#ax1.gca().patch.set_color('xkcd:tan') #'xkcd:olive' 'xkcd:khaki' 'xkcd:green'
#ax1.contour(z=tau_scene,cmap='Reds',corner_mask=corner_mask)
#a = ax1.colorbar()
#a.set_label('Cloud Optical Thickness')
#
#ax2  = axs[1]
#ax2.gca().patch.set_color('xkcd:tan') #'xkcd:olive' 'xkcd:khaki' 'xkcd:green'
#ax2.contour(z=reff_scene,cmap='Blues',corner_mask=corner_mask)
#b = ax2.colorbar()
#b.set_label('Effective Radius (um)')