# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 19:18:05 2019

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

#phase_dict=water_dict
phase_dict=ice_dict

solar_dict = read_solar(solar_path)
hysics_dict = read_HySICS(phase_dict) #AAAAAAAAAAAAAA  
hysics_dict = shift_hysics(hysics_dict)
hysics_dict = calc_HySICS_reflectance(hysics_dict,solar_dict)
LRT_dict = read_all_LRT(phase_dict) #AAAAAAAAAAAAAA    

plot_all_LRT(LRT_dict,hysics_dict)    


#For Many Retrieval Algorithms
nums=np.arange(0,20)
dict_of_analysis_dicts={}
for n in nums:
    wl_num = (n+1)*5
    analysis_dict = disagreement_all_LRT(hysics_dict, wl_num, phase_dict)
    dict_of_analysis_dicts['analysis_dict_'+str(wl_num)] = analysis_dict

for k in dict_of_analysis_dicts.keys():
    d = dict_of_analysis_dicts[k]
    chi_sq_contour(d)
    plt.xlim(5,60) #30
    plt.ylim(10,170) #plt.ylim(2,10)

plot_algorithm_comparison(dict_of_analysis_dicts) 

#LRT_best_dict = plot_best_fits(hysics_dict, LRT_dict, dict_of_analysis_dicts)
#plt.ylim(10**-4,1)
#plt.xlim(740,780)

#Next steps,
#Pick best retrieval and apply to every pixel, make plot of reff and COT
#Do I need to repeat LRT ice output for baum verbose
#More reff and tau's around retrieved value for ice cloud modeling
