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

data_cube = np.load('data_cube_water.npy')
rgb = np.load(phase_dict['rgb_path'])

xx = len(data_cube[:,0,0])
yy = len(data_cube[0,:,0])

tau_scene = np.zeros((xx,yy), dtype = float)
reff_scene = np.zeros((xx,yy), dtype = float)
mean_scene = np.zeros((xx,yy), dtype = float)
mask = np.zeros((xx,yy), dtype = float)

for y in tqdm(np.arange(0,yy)):
    for x in np.arange(0,xx):
        mean = np.mean(rgb[x,y])
        if (mean<0.75): #No Clouds
                tau_scene[x,y] = 'NaN'
                reff_scene[x,y] = 'NaN'
        elif:
                pixel = data_cube[x,y,:]
                hysics_data = pixel[0:-27] if phase_dict['phase'] == 'Liquid Water' else pixel[35:-27]
                hysics_dict = {'wl':hysics_wl,'data':pixel,'phase':phase_dict['phase']}
                
                hysics_dict= calc_HySICS_reflectance(hysics_dict,solar_dict)
                analysis_dict = disagreement_all_LRT(hysics_dict, wl_num, phase_dict) 
                
                tau_scene[x,y] = analysis_dict['best_COT']
                reff_scene[x,y] = analysis_dict['best_reff']


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