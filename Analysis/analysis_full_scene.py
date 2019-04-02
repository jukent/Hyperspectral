import numpy as np
import matplotlib.pyplot as plt
from help_analysis_read import read_all_LRT, read_solar
from hel_analysis import calc_HySICS_reflectance, disagreement_all_LRT
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

for x in np.arange(0,xx):
    for y in np.arange(0,yy):
        hysics_data = [x,y,:]
        hysics_dict = {'wl':wl_list,'data':hysics_data,'phase':phase_dict['phase']}
        rgb_pixel = [x,y,:]

        hysics_dict = calc_HySICS_reflectance(hysics_dict,solar_dict)
        hysics_dict = mask_clouds(hysics_dict,rgb_pixel)

        analysis_dict = disagreement_all_LRT(hysics_dict, wl_num, phase_dict) 
        tau_scene[x,y] = analysis_dict['best_COT']
        reff_scene[x,y] = analysis_dict['best_reff']

#Cloud Mask
def mask_clouds(hysics_dict,rgb_pixel):
    data = hysics_dict['data']
    if (np.abs(rgb[0]/rgb[1] - rgb[2]/rgb[1]) > 0.2):
        [i='NaN' for i in data]
    hysics_dict['data']=data
    return hysics_dict;

#Plot tau_scene and reff_scene -- Needs Testing
plt.scatter(data_cube[:,0,0],data_cube[0,:,0],c=tau_scene)
plt.scatter(data_cube[:,0,0],data_cube[0,:,0],c=reff_scene)