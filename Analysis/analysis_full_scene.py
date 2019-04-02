import numpy as np
import matplotlib.pyplot as plt
from help_analysis_read import read_all_LRT, read_solar
from hel_analysis import calc_HySICS_reflectance, disagreement_all_LRT, mask_ground
from help_analysis_plot import plot_all_LRT, chi_sq_contour, plot_algorithm_comparison, plot_best_fits
import generate_rgb_cube

solar_path = '../../../CU Boulder/LASP/LASP py/Solar Data/Aug14SolarRef.dat'
HySICS_wl_path = '../../../CU Boulder/LASP/LASP py/HySICS Data/WLHysics.sav'
HySICS_LW_data_path = '../../../CU Boulder/LASP/LASP py/Desert_vegetation_clouds/'
HySICS_IC_data_path = '../../../CU Boulder/LASP/LASP py/Thick_clouds1'
LRT_LW_path = '../../../CU Boulder/LASP/LASP py/WC files/Outputs/'
verbose_LW_path = '../../../CU Boulder/LASP/LASP py/WC files/Verbose/'
LRT_IC_path = '../../../CU Boulder/LASP/LASP py/ICfiles/Repeat for Baum/Outputs/ghm/'#solid-column/'#rough-aggregate/'
verbose_IC_path = '../../../CU Boulder/LASP/LASP py/ICfiles/Verbose/ghm/'#solid-column/'#rough-aggregate/'

water_dict = {'solar_path':solar_path,'HySICS_wl_path':HySICS_wl_path,'HySICS_data_path':HySICS_LW_data_path, \
              'LRT_path':LRT_LW_path,'LRT_verbose_path':verbose_LW_path, 'phase':'Liquid Water'}

ice_dict = {'solar_path':solar_path,'HySICS_wl_path':HySICS_wl_path,'HySICS_data_path':HySICS_IC_data_path, \
              'LRT_path':LRT_IC_path,'LRT_verbose_path':verbose_IC_path, 'phase':'Ice'}

solar_dict = read_solar(solar_path)

#-------------------------
wl_num = 15
phase_dict = water_dict
#phase_dict = ice_dict
#-------------------------

#Data_cube HySICS values are shifted by 4 nm
(data_cube, wl_list, rgb) = generate_rgb_cube(phase_dict['HySICS_data_path'])

xx = len(data_cube[:,0,0])
yy = len(data_cube[0,:,0])
tau_scene = zeros(xx,yy)
reff_scene = zeros(xx,yy)
mask = zeros(xx,yy)

for x in np.arange(0,xx):
        for y in np.arange(0,yy):
                mask_ground[x,y]  = mask_ground(rgb[x,y,:])

                if (mask_ground[x,y] == True): #Ground -- No Clouds
                        tau_scene[x,y] = 'NaN'
                        reff_scene[x,y] = 'NaN'
                        
                elif (mask_ground[x,y] == False): #No Ground -- Yes Clouds
                        hysics_data = [x,y,:]
                        hysics_dict = {'wl':wl_list,'data':hysics_data,'phase':phase_dict['phase']}
                        hysics_dict= calc_HySICS_reflectance(hysics_dict,solar_dict)
                        
                        analysis_dict = disagreement_all_LRT(hysics_dict, wl_num, phase_dict) 
                        tau_scene[x,y] = analysis_dict['best_COT']
                        reff_scene[x,y] = analysis_dict['best_reff']


#Check Mask
fig, axs = plt.subplots(1,2,sharex=True)
ax1 = axs[0]
ax1.imshow(rgb)

rgb_masked = rgb*mask_ground
ax2 = axs[1]
ax2.imshow(rgb_masked)


#Plot tau_scene and reff_scene -- Needs Testing
fig, axs = plt.subplots(1,2,sharex=True)

ax1  = axs[0]
ax1.gca().patch.set_color('xkcd:tan') #'xkcd:olive' 'xkcd:khaki' 'xkcd:green'
ax1.contour(z=tau_scene,cmap='Reds',corner_mask=corner_mask)
a = ax1.colorbar()
a.set_label('Cloud Optical Thickness')

ax2  = axs[1]
ax2.gca().patch.set_color('xkcd:tan') #'xkcd:olive' 'xkcd:khaki' 'xkcd:green'
ax2.contour(z=reff_scene,cmap='Blues',corner_mask=corner_mask)
b = ax2.colorbar()
b.set_label('Effective Radius (um)')