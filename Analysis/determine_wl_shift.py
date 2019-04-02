# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 17:24:16 2019

@author: julia
"""
import numpy as np
from help_read import read_HySICS, read_all_LRT, read_solar
from helpers import calc_HySICS_reflectance, find_retrieval_wavelengths, disagreement_all_LRT
from help_plot import plot_best_fits


solar_path = '../../../CU Boulder/LASP/LASP py/Solar Data/Aug14SolarRef.dat'
HySICS_wl_path = '../../../CU Boulder/LASP/LASP py/HySICS Data/WLHysics.sav'
HySICS_LW_data_path = '../../../CU Boulder/LASP/LASP py/WC files/thincloud.npy'
LRT_LW_path = '../../../CU Boulder/LASP/LASP py/WC files/Outputs/'
verbose_LW_path = '../../../CU Boulder/LASP/LASP py/WC files/Verbose/'

water_dict = {'solar_path':solar_path,'HySICS_wl_path':HySICS_wl_path,'HySICS_data_path':HySICS_LW_data_path, \
              'LRT_path':LRT_LW_path,'LRT_verbose_path':verbose_LW_path, 'phase':'Liquid Water'}

solar_dict = read_solar(solar_path)
hysics_dict = read_HySICS(water_dict) 
hysics_dict = calc_HySICS_reflectance(hysics_dict,solar_dict)
LRT_dict = read_all_LRT(water_dict)

wl_num=15 #must be a factor of 5
wl_list = find_retrieval_wavelengths(wl_num) 

analysis_dict = disagreement_all_LRT(hysics_dict, wl_num, water_dict) #AAAAAAAAAAA
dict_of_analysis_dicts={'analysis_dict':analysis_dict}
LRT_best_dict = plot_best_fits(hysics_dict, LRT_dict, dict_of_analysis_dicts)


LRT_idx_min = np.abs([x-750 for x in LRT_dict['wl']]).argmin()
LRT_idx_max = np.abs([x-780 for x in LRT_dict['wl']]).argmin()
LRT_val_filt = LRT_best_dict['analysis_dict'][LRT_idx_min:LRT_idx_max]
LRT_wl_filt = LRT_dict['wl'][LRT_idx_min:LRT_idx_max]

HySICS_idx_min = np.abs([x-750 for x in hysics_dict['wl']]).argmin()
HySICS_idx_max = np.abs([x-780 for x in hysics_dict['wl']]).argmin()
HySICS_val_filt = hysics_dict['data'][HySICS_idx_max:HySICS_idx_min]
HySICS_wl_filt = hysics_dict['wl'][HySICS_idx_max:HySICS_idx_min]

idx_LRT_O2 = np.abs(LRT_val_filt).argmin()
idx_HySICS_O2 = np.abs(HySICS_val_filt).argmin()

wl_LRT_O2 = LRT_wl_filt[idx_LRT_O2]
wl_HySICS_O2 = HySICS_wl_filt[idx_HySICS_O2]

shift = wl_LRT_O2 - wl_HySICS_O2 #5.77

def shift_hysics(hysics_dict):
    wl = [x+5.77 for x in hysics_dict['wl']]
    hysics_dict['wl']=wl
    return hysics_dict;