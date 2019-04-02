# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 02:28:23 2019

@author: julia
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def plot_all_LRT(LRT_dict,hysics_dict):
    plt.figure()
    plt.xlim([450,2200])
    plt.xlabel('Wavlength (nm)',fontsize=18)
    plt.ylabel('Radiance (W/m^2/nm/sr)',fontsize=18)
    tt = 'Thin Water Cloud Radiance' if hysics_dict['phase']=='Liquid Water' else 'Thick Ice Cloud Radiance'
    plt.title(tt,fontsize=25)
    label_size = 20
    plt.rcParams['xtick.labelsize'] = label_size
    plt.rcParams['ytick.labelsize'] = label_size
    for d in LRT_dict['data']:
        plt.semilogy(LRT_dict['wl'],d,'--',c='lightgray')
    plt.semilogy(hysics_dict['wl'],hysics_dict['data'],'black')
    return;


def chi_sq_contour(analysis_dict):           
    plt.figure()
    
    cnorm = mpl.colors.Normalize(vmin=0.,vmax=20.)#40.)
    
    x = analysis_dict['r_eff']
    y = analysis_dict['COT']
    z = analysis_dict['difference']
    b = plt.tricontour(x,y,z,70,norm=cnorm)#70)
    
    plt.xlabel('Cloud Particle Effective Radius (r_eff) (um)',fontsize=18)
    plt.ylabel('Cloud Optical Thickness (COT)',fontsize=18)
    plt.title(analysis_dict['phase']+' Cloud {} Wavelength Retrieval Values'.format(str(analysis_dict['num_retrieval_wavelengths'])),fontsize=20)
    plt.clabel(b, inline=1, fontsize=15)
    
    reff_best = analysis_dict['best_reff']
    COT_best = analysis_dict['best_COT'] 
    label = 'COT: '+str(round(COT_best,2))+', r_eff: '+str(reff_best)+'um'
    bbox = dict(boxstyle="round", fc="0.8",alpha=1)
    arrow = dict(facecolor='black', shrink=0.05)
    plt.annotate(label,[reff_best,COT_best],xytext=(reff_best+2,COT_best+0.5),\
        bbox=bbox,fontsize=15,arrowprops=arrow)
    plt.scatter(reff_best,COT_best,marker="X",c='k',s=70)
    return;
   
    
def plot_algorithm_comparison(dict_of_analysis_dicts):
    best_COTs = []
    best_reffs = []
    wl_nums = []
    for k in dict_of_analysis_dicts.keys():
        d = dict_of_analysis_dicts[k]
        
        best_COTs.append(d['best_COT'])
        best_reffs.append(d['best_reff'])
        wl_nums.append(d['num_retrieval_wavelengths'])
    
    fig, ax1 = plt.subplots()
    
    color = 'tab:red'
    ax1.set_xlabel('Number of Retrieval Wavelengths', fontsize=18)
    ax1.set_ylabel('Cloud Optical Thickness (COT)', color=color, fontsize=18)
    ax1.plot(wl_nums, best_COTs, color=color)
    ax1.scatter(wl_nums, best_COTs,c=color)
    ax1.tick_params(axis='y', labelcolor=color)
    
    ax2 = ax1.twinx()
    
    color = 'tab:blue'
    ax2.set_ylabel('Cloud Particle Effective Radius (r_eff) (um)', color=color, fontsize=18)
    ax2.plot(wl_nums, best_reffs, color=color)
    ax2.scatter(wl_nums, best_reffs, c=color)
    ax2.tick_params(axis='y', labelcolor=color)
    
    fig.tight_layout()
    plt.title('Cloud Microphysical Parameters from Different '+ \
        'Retrieval Algorithms for '+d['phase']+' Cloud',fontsize=25)
    plt.show()
    return;
    

def plot_best_fits(hysics_dict, LRT_dict, dict_of_analysis_dicts):            
    plt.figure()
    plt.xlim([440,2210])
    plt.ylim(top=0.3)
    plt.xlabel('Wavlength (nm)',fontsize=18)
    plt.ylabel('Radiance (W/m^2/nm/sr)',fontsize=18)
    plt.title(hysics_dict['phase']+' Cloud -  Best Fit',fontsize=25)
#    plt.title('Demonstration of HySICS Shift (unshifted)',fontsize=25)
            
    plt.plot(hysics_dict['wl'],hysics_dict['data'],'black',label='HySICS')
    
    LRT_best_dict = {}
    for k in dict_of_analysis_dicts.keys():
        d = dict_of_analysis_dicts[k]
        
        LRT_best = LRT_dict['data'][d['best_index']]
        LRT_best_dict[k] = LRT_best
        
        label = 'LibRadtran'
#        label ='COT: '+str(round(d['best_COT'],2))+ \
#            ' r_eff: '+str(d['best_reff'])+'um -- '+ \
#                str(d['num_retrieval_wavelengths'])+' Wavelengths Algorithm: '+ \
#                    '%.2f'%d['best_difference']+'%'        
        plt.semilogy(LRT_dict['wl'],LRT_best,label=label)        
        [plt.axvline(w,linestyle='--',c='grey') for w in np.ravel(d['retrieval_wavelengths'])]
#        [plt.axvline(w,linestyle='--',c='grey') for w in [762]]
    
    plt.legend(fontsize=20)
    plt.ylim([0,0.15])
    label_size = 20
    plt.rcParams['xtick.labelsize'] = label_size
    plt.rcParams['ytick.labelsize'] = label_size
    return (LRT_best_dict);