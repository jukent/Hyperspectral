# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 19:18:05 2019

@author: julia
"""

import numpy as np
from help_read import read_HySICS, read_all_LRT, read_solar
from helpers import calc_HySICS_reflectance, find_retrieval_wavelengths, disagreement_all_LRT
from help_plot import plot_all_LRT, chi_sq_contour, plot_algorithm_comparison, plot_best_fits


solar_path = '../../../CU Boulder/LASP/LASP py/Solar Data/Aug14SolarRef.dat'

HySICS_wl_path = '../../../CU Boulder/LASP/LASP py/HySICS Data/WLHysics.sav'
HySICS_LW_data_path = '../../../CU Boulder/LASP/LASP py/WC files/thincloud.npy'
HySICS_IC_data_path = '../../../CU Boulder/LASP/LASP py/ICfiles/thickcloud.npy'

LRT_LW_path = '../../../CU Boulder/LASP/LASP py/WC files/Outputs/'
verbose_LW_path = '../../../CU Boulder/LASP/LASP py/WC files/Verbose/'
LRT_IC_path = '../../../CU Boulder/LASP/LASP py/ICfiles/Repeat for Baum/Outputs/ghm/'
verbose_IC_path = '../../../CU Boulder/LASP/LASP py/ICfiles/Verbose/ghm/'

water_dict = {'solar_path':solar_path,'HySICS_wl_path':HySICS_wl_path,'HySICS_data_path':HySICS_LW_data_path, \
              'LRT_path':LRT_LW_path,'LRT_verbose_path':verbose_LW_path, 'phase':'Liquid Water'}

ice_dict = {'solar_path':solar_path,'HySICS_wl_path':HySICS_wl_path,'HySICS_data_path':HySICS_IC_data_path, \
              'LRT_path':LRT_IC_path,'LRT_verbose_path':verbose_IC_path, 'phase':'Ice'}


hysics_dict = read_HySICS(ice_dict)
LRT_dict = read_all_LRT(ice_dict)

solar_dict = read_solar(solar_path)
hysics_dict, solar_dict = calc_HySICS_reflectance(hysics_dict,solar_dict)


        
#Plot all LRT and HySICS        
plot_all_LRT(LRT_dict,hysics_dict)    
    
#For One Retrieval Algorithm
wl_num=15 #must be a factor of 5
wl_list = find_retrieval_wavelengths(wl_num)      
analysis_dict = disagreement_all_LRT(hysics_dict, solar_dict, wl_num, wl_list, ice_dict)
chi_sq_contour(analysis_dict)


##For Many Retrieval Algorithms
#nums=np.arange(0,3)
#dict_of_analysis_dicts={}
#for n in nums:
#    wl_num = (n+1)*5
#    wl_list = find_retrieval_wavelengths(wl_num) 
#    analysis_dict = disagreement_all_LRT(hysics_dict, solar_dict, wl_num, wl_list, LRT_LW_path, verbose_LW_path)
#    dict_of_analysis_dicts['analysis_dict_'+str(wl_num)] = analysis_dict
#
##Plot chi_squared residuals for all number algorithms
#for k in dict_of_analysis_dicts.keys():
#    d = dict_of_analysis_dicts[k]
#    chi_sq_contour(d)
#
##Plot tau vs wl_num and reff vs wl_num
#plot_algorithm_comparison(dict_of_analysis_dicts) 
#
##Plot Best Fits
#plot_best_fits(hysics_dict, LRT_dict, dict_of_analysis_dicts)



#Next steps,
#HySICs shift
#Pick best retrieval and apply to every pixel, make plot of reff and COT

#Send just 15 channel csv and plot for water and ice to Peter