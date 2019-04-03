# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 11:56:10 2019

@author: julia
"""

import numpy as np
import matplotlib.pyplot as plt
from help_read import read_HySICS, read_all_LRT, read_solar
from helpers import calc_HySICS_reflectance, disagreement_all_LRT
from help_plot import plot_all_LRT, chi_sq_contour, plot_algorithm_comparison, plot_best_fits

def shift_hysics(hysics_dict):
    wl = [x+4 for x in hysics_dict['wl']]
    hysics_dict['wl']=wl
    return hysics_dict;    

solar_path = '../../../CU Boulder/LASP/LASP py/Solar Data/Aug14SolarRef.dat'
HySICS_wl_path = '../../../CU Boulder/LASP/LASP py/HySICS Data/WLHysics.sav'
HySICS_LW_data_path = '../../../CU Boulder/LASP/LASP py/WC files/thincloud.npy'
HySICS_IC_data_path = '../../../CU Boulder/LASP/LASP py/ICfiles/thickcloud.npy'
LRT_LW_path = '../../../CU Boulder/LASP/LASP py/WC files/Outputs/'
verbose_LW_path = '../../../CU Boulder/LASP/LASP py/WC files/Verbose/'
LRT_IC_path = '../../../CU Boulder/LASP/LASP py/ICfiles/Repeat for Baum/Outputs/ghm/'#solid-column/'#rough-aggregate/'
verbose_IC_path = '../../../CU Boulder/LASP/LASP py/ICfiles/Verbose/ghm/'#solid-column/'#rough-aggregate/'

water_dict = {'solar_path':solar_path,'HySICS_wl_path':HySICS_wl_path,'HySICS_data_path':HySICS_LW_data_path, \
              'LRT_path':LRT_LW_path,'LRT_verbose_path':verbose_LW_path, 'phase':'Liquid Water'}

ice_dict = {'solar_path':solar_path,'HySICS_wl_path':HySICS_wl_path,'HySICS_data_path':HySICS_IC_data_path, \
              'LRT_path':LRT_IC_path,'LRT_verbose_path':verbose_IC_path, 'phase':'Ice'}

solar_dict = read_solar(solar_path)

#Water Unshifted
phase_dict=water_dict

hysics_dict = read_HySICS(phase_dict) #AAAAAAAAAAAAAA  
hysics_dict = calc_HySICS_reflectance(hysics_dict,solar_dict)
LRT_dict = read_all_LRT(phase_dict) #AAAAAAAAAAAAAA    

nums=np.arange(2,3)
dict_of_analysis_dicts={}
for n in nums:
    wl_num = (n+1)*5
    analysis_dict = disagreement_all_LRT(hysics_dict, wl_num, phase_dict)
    dict_of_analysis_dicts['analysis_dict_'+str(wl_num)] = analysis_dict

LRT_best_dict = plot_best_fits(hysics_dict, LRT_dict, dict_of_analysis_dicts)
plt.ylim(10**-4,1)

HySICS_water_unshifted_dict = {'Radiances':hysics_dict['data'],'Wavelengths':hysics_dict['wl']}
LRT_water_unshifted_dict = {'Radiances':LRT_best_dict['analysis_dict_15'],'Wavelengths':LRT_dict['wl']}

f = open('../../HySICS_water_unshifted.txt','w')
f.write('HySICS Thin Cloud - Unshifted' + repr(HySICS_water_unshifted_dict))
f.close()
f = open('../../LRT_water_unshifted.txt','w')
f.write('LRT Thin Cloud - Unshifted' + repr(LRT_water_unshifted_dict))
f.close()

#Water Shifted
hysics_dict = read_HySICS(phase_dict) #AAAAAAAAAAAAAA  
hysics_dict = shift_hysics(hysics_dict)
hysics_dict = calc_HySICS_reflectance(hysics_dict,solar_dict)
LRT_dict = read_all_LRT(phase_dict) #AAAAAAAAAAAAAA    

nums=np.arange(2,3)
dict_of_analysis_dicts={}
for n in nums:
    wl_num = (n+1)*5
    analysis_dict = disagreement_all_LRT(hysics_dict, wl_num, phase_dict)
    dict_of_analysis_dicts['analysis_dict_'+str(wl_num)] = analysis_dict

LRT_best_dict = plot_best_fits(hysics_dict, LRT_dict, dict_of_analysis_dicts)
plt.ylim(10**-4,1)

HySICS_water_shifted_dict = {'Radiances':hysics_dict['data'],'Wavelengths':hysics_dict['wl']}
LRT_water_shifted_dict = {'Radiances':LRT_best_dict['analysis_dict_15'],'Wavelengths':LRT_dict['wl']}

f = open('../../HySICS_water_shifted.txt','w')
f.write('HySICS Thin Cloud - Shifted' + repr(HySICS_water_shifted_dict))
f.close()
f = open('../../LRT_water_shifted.txt','w')
f.write('LRT Thin Cloud - Shifted' + repr(LRT_water_shifted_dict))
f.close()


#Ice Unshifted
phase_dict=ice_dict

hysics_dict = read_HySICS(phase_dict) #AAAAAAAAAAAAAA  
hysics_dict = calc_HySICS_reflectance(hysics_dict,solar_dict)
LRT_dict = read_all_LRT(phase_dict) #AAAAAAAAAAAAAA    

nums=np.arange(2,3)
dict_of_analysis_dicts={}
for n in nums:
    wl_num = (n+1)*5
    analysis_dict = disagreement_all_LRT(hysics_dict, wl_num, phase_dict)
    dict_of_analysis_dicts['analysis_dict_'+str(wl_num)] = analysis_dict

LRT_best_dict = plot_best_fits(hysics_dict, LRT_dict, dict_of_analysis_dicts)
plt.ylim(10**-4,1)

HySICS_ice_unshifted_dict = {'Radiances':hysics_dict['data'],'Wavelengths':hysics_dict['wl']}
LRT_ice_unshifted_dict = {'Radiances':LRT_best_dict['analysis_dict_15'],'Wavelengths':LRT_dict['wl']}

f = open('../../HySICS_ice_unshifted.txt','w')
f.write('HySICS Thick Cloud - Unshifted' + repr(HySICS_ice_unshifted_dict))
f.close()
f = open('../../LRT_ice_unshifted.txt','w')
f.write('LRT Thick Cloud - Unshifted' + repr(LRT_ice_unshifted_dict))
f.close()

#Ice Shifted
hysics_dict = read_HySICS(phase_dict) #AAAAAAAAAAAAAA  
hysics_dict = shift_hysics(hysics_dict)
hysics_dict = calc_HySICS_reflectance(hysics_dict,solar_dict)
LRT_dict = read_all_LRT(phase_dict) #AAAAAAAAAAAAAA    

nums=np.arange(2,3)
dict_of_analysis_dicts={}
for n in nums:
    wl_num = (n+1)*5
    analysis_dict = disagreement_all_LRT(hysics_dict, wl_num, phase_dict)
    dict_of_analysis_dicts['analysis_dict_'+str(wl_num)] = analysis_dict

LRT_best_dict = plot_best_fits(hysics_dict, LRT_dict, dict_of_analysis_dicts)
plt.ylim(10**-4,1)

HySICS_ice_shifted_dict = {'Radiances':hysics_dict['data'],'Wavelengths':hysics_dict['wl']}
LRT_ice_shifted_dict = {'Radiances':LRT_best_dict['analysis_dict_15'],'Wavelengths':LRT_dict['wl']}

f = open('../../HySICS_ice_shifted.txt','w')
f.write('HySICS Thick Cloud - Shifted' + repr(HySICS_ice_shifted_dict))
f.close()
f = open('../../LRT_ice_shifted.txt','w')
f.write('LRT Thick Cloud - Shifted' + repr(LRT_ice_shifted_dict))
f.close()