# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 23:15:49 2019

@author: julia
"""
from scipy.interpolate import interp1d
import numpy as np
import os
import help_analysis_read as rd

def calc_HySICS_reflectance(hysics_dict,solar_dict):    
    f=interp1d(solar_dict['wl'],solar_dict['flux'])
    flux_sol_new = f(hysics_dict['wl'])
    
    refl_HySICS= [np.pi*hysics_dict['data'][w]/flux_sol_new[w] \
        for w in np.arange(0,len(hysics_dict['wl']))]
    
    hysics_dict['refl'] = refl_HySICS
    hysics_dict['solar_flux'] = flux_sol_new
    return hysics_dict;


def calc_LRT_reflectance(hysics_dict,file):   
    (x,y) = rd.read_LRT(file,hysics_dict['phase'])
    f = interp1d(x,y)
    xnew = f(hysics_dict['wl'][1:])
    refl_LRT = [np.pi*xnew[w]/hysics_dict['solar_flux'][w] \
        for w in np.arange(0,len(hysics_dict['wl'][1:]))]
    return refl_LRT;


def find_retrieval_wavelengths(num_wl): 
    if num_wl%5 !=0: print('Number of wavelengths in algorithm must be a factor of 5.')
    wl_lims = [[600,750], [980,1050], [1220,1320], [1500,1750], [2100,2200]]
    wl_list = []
    num = num_wl/5
    for w in np.arange(0,5):
        wl = np.linspace(wl_lims[w][0],wl_lims[w][1],num)
        wl_list.append(wl)
    return wl_list;
    if (num_wl ==5):
        wl_list=[[750], [1000], [1200], [1660], [2200]]

def disagreement_algorithm(hysics_dict, wl_num, wl_list, path, file, verbose_path):  
    refl_LRT = calc_LRT_reflectance(hysics_dict,path+file)
    idx_1 = (np.abs(hysics_dict['wl'] - wl_list[0][0])).argmin()
    
    sigma_all = []
    for i in np.arange(0,5):
        k=i+1
        
        sigma_k = []      
        for j in np.arange(0,len(wl_list[i])):
            idx_wl = (np.abs(hysics_dict['wl']-wl_list[i][j])).argmin()
            
            reflectance_term = (5-k)**2*(refl_LRT[idx_wl]-hysics_dict['refl'][idx_wl])**2
            ratio_term = (k-1)**2*(refl_LRT[idx_wl]/refl_LRT[idx_1]- \
                hysics_dict['refl'][idx_wl]/hysics_dict['refl'][idx_1])**2
            
            sigma_k_j = reflectance_term + ratio_term
            sigma_k.append(sigma_k_j)
        sigma_all.append(sigma_k)
        
    sigma = np.sum(sigma_all)
    chi_sq = sigma/(wl_num*12)
    
    diff = np.sqrt(chi_sq)*100
    
    fs = file.split('_')
    reff = fs[3][1:] if hysics_dict['phase']=='Liquid Water' else fs[2][1:]
    ff = file[0:-4] if hysics_dict['phase']=='Liquid Water' else file[0:-9]
    verbosefile = verbose_path+ff+'_550nm_verbose.txt'
    tau = rd.read_verbose(verbosefile,hysics_dict['phase'])
    
    diff_dict = {'phase':hysics_dict['phase'],'r_eff':[reff],'COT':tau,'difference':diff,\
        'wavelengths':hysics_dict['wl'],'num_retrieval_wavelengths':wl_num,'retrieval_wavelengths':wl_list}
    return diff_dict;


def disagreement_all_LRT(hysics_dict, wl_num, phase_dict):
    wl_list = find_retrieval_wavelengths(wl_num) 
    
    reff_list = [] ; COT_list = [] ; diff_list = []
    for f in os.listdir(phase_dict['LRT_path']):
        if f.endswith(".dat"):
            diff_dict = disagreement_algorithm(hysics_dict, wl_num, wl_list, phase_dict['LRT_path'], f, phase_dict['LRT_verbose_path'])
            reff_list.append(diff_dict['r_eff'])
            COT_list.append(diff_dict['COT'])   
            diff_list.append(diff_dict['difference'])
    val_diff, idx_diff = min((val, idx) for (idx, val) in enumerate(diff_list))

    reff_float = [float(i) for sublist in reff_list for i in sublist]
    COT_float = [float(i) for i in COT_list]
    reff_best = float(reff_float[idx_diff])
    COT_best = float(COT_float[idx_diff])  
    
    analysis_dict = {'phase':hysics_dict['phase'],'r_eff':reff_float,'COT':COT_float,'difference':diff_list, \
        'wavelengths':diff_dict['wavelengths'],'num_retrieval_wavelengths':wl_num,'retrieval_wavelengths':wl_list, \
            'best_index':idx_diff,'best_COT':COT_best,'best_reff':reff_best,'best_difference':val_diff}
    return(analysis_dict);

    
def shift_hysics(hysics_dict):
    wl = [x+4 for x in hysics_dict['wl']]
    hysics_dict['wl']=wl
    return hysics_dict;    
        
