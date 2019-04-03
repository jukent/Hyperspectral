import dask.dataframe as dd
import dask.array as da
import os
import xarray as xr
import csv
import numpy as np
import pandas as pd
import dask

solar_path = 'Aug14SolarRef.dat'
HySICS_wl_path = 'HySICS files/WLHysics.sav'

HySICS_LW_data_path = 'Desert_vegetation_clouds/'
LW_datacube_path = 'data_cube_water.npy'
LW_rgb_path = 'rgb_water.npy'
LRT_LW_path = 'WC files/'
verbose_LW_path = 'WC files/verbose/'

HySICS_IC_data_path = 'Thick_clouds1/'
IC_datacube_path = 'data_cube_ice.npy'
IC_rgb_path = 'rgb_ice.npy'
LRT_IC_path = 'IC files/ghm/'
verbose_IC_path = 'IC files/ghm/verbose_ghm/'

water_dict = {'solar_path':solar_path,'HySICS_wl_path':HySICS_wl_path,'HySICS_data_path':HySICS_LW_data_path, \
              'LRT_path':LRT_LW_path,'LRT_verbose_path':verbose_LW_path, 'phase':'Liquid Water', \
             'data_cube_path':LW_datacube_path,'rgb_path':LW_rgb_path}
ice_dict = {'solar_path':solar_path,'HySICS_wl_path':HySICS_wl_path,'HySICS_data_path':HySICS_IC_data_path, \
              'LRT_path':LRT_IC_path,'LRT_verbose_path':verbose_IC_path, 'phase':'Ice',\
           'data_cube_path':IC_datacube_path,'rgb_path':IC_rgb_path}

phase_dict=water_dict           

def read_HySICS_wl(phase_dict):
    hysics_wl = sio.readsav(phase_dict['HySICS_wl_path'])
    hysics_wl = hysics_wl['wlsample']
    hysics_wl = hysics_wl[0:-27] if phase_dict['phase']=='Liquid Water' else hysics_wl[35:-27]
    hysics_wl = [x+4 for x in hysics_wl]
    return hysics_wl;
wl_list = read_HySICS_wl

@dask.delayed
def read_df(file):
    array = np.genfromtxt(file,delimiter=',')
    return array;

def generate_datacube(phase_dict,wl_list)
    path = phase_dict['HySICS_data_path']
    files = [os.path.join(path,f) for f in os.listdir(path) if f.startswith('img')]
    files.sort()

    dfs = [read_df(file) for file in files]

    sample = dfs[0].compute()
    das = [da.from_delayed(item, shape=sample.shape, dtype=sample.dtype) for item in dfs]
    array = da.stack(das)

    datacube = xr.DataArray(array, dims=['x', 'y', 'wavelength'], name='datacube',coords=['wavelength':wl_list])
    return datacube;

datacube_dvc = generate_datacube(phase_dict,wl_list)