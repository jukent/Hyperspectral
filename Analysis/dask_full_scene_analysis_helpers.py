import csv
import numpy as np
import xarray as xr
from scipy.interpolate import interp1d
import os


def read_verbose(verbose_file,phase):
    data = [x for x in open(verbose_file,'r')]
    a = data.index('Using new intensity correction, with phase functions\n')
    n=4 if phase=='Liquid Water' else 5
    tau = data[a:][67].split('|')[n].split()[0]
    return(tau);


def read_LRT(file,phase):
    fs = file.split('_')
    reff = fs[3][1:] if phase =='Liquid Water' else fs[2][1:]
    
    fv = file.split('/')
    ff = fv[1][0:-4] if phase =='Liquid Water' else fv[1][0:-9]
    verbose_file = fv[0]+'/verbose/'+ff+'_550nm_verbose.txt'
    tau = read_verbose(verbose_file,phase)
    
    data = [x for x in csv.reader(open(file,'r'),delimiter='\t')]  
    data = data[5:]
    wl=[]           #Wavelength (nm)
    dirL=[]         #Direct Radiance (W/m^2/nm/sr)
    
    loops = int(len(data)/3)
    if phase =='Liquid Water': 
        dirL = [np.float(data[3*n+1][0][1:9]) for n in range(loops)]
        wl = [np.float(data[3*n+1][0][1:9]) for n in range(loops)]
    else:
        dirL = [np.float(data[3*n+2][0][9:24]) for n in range(loops)]
        wl = [np.float(data[3*n][0][1:9]) for n in range(loops)]

    dirLW = [i*10**-3 for i in dirL] #convert to W/m^2/nm/sr
    
    lrt_radiances = xr.DataArray(dirLW, dims=['wavelength'], name='LRT', coords={'wavelength':wl})
    lrt_radiances.attrs['COT'] = tau
    lrt_radiances.attrs['r_eff'] = reff
    return lrt_radiances;


def read_solar(file):
    d = np.loadtxt(file, delimiter="\t")
    solar_wl=d[:,0]
    solar_flux=d[:,1]
    solar_flux = np.array([i/1000 for i in solar_flux]) 
    solar = xr.DataArray(solar_flux,dims=['wavelength'], name='solar flux', coords={'wavelength':solar_wl})
    return solar;    


def filter_LRT_wl(lrt_rad):
    wl_min = 450 ; wl_max = 2300
    lrt_radf = lrt_rad.where(lrt_rad.wavelength>wl_min, drop=True)
    lrt_radf = lrt_radf.where(lrt_radf.wavelength<wl_max, drop=True)
    return lrt_radf;


def get_sample_lrt_radf(path,phase):
    files = [os.path.join(path,f) for f in os.listdir(path)]
    file = files[0]
    lrt_rad = read_LRT(file,phase)
    lrt_radf_samp = filter_LRT_wl(lrt_rad)
    return (files, lrt_radf_samp);


def interpolate_solar(lrt_radf,solar):
    f = interp1d(solar.wavelength,solar.values,fill_value='NaN')
    solar_new = f(lrt_radf.wavelength)
    
    solar_interp = xr.DataArray(solar_new,dims=['wavelength'], name='interpolated solar flux', coords={'wavelength':lrt_radf.wavelength})
    return solar_interp;


def calc_LRT_reflectance(lrt_radf,solar_interp,phase):
    mu = np.cosd(38.5) if phase == 'Liquid Water' else np.cosd(56)
    refl_lrt = [np.pi/mu*lrt_radf.values[w]/solar_interp.values[w] for w in np.arange(0,len(lrt_radf.wavelength))]
                                                                        
    refl_lrt = xr.DataArray(refl_lrt,dims=['wavelength'], name='lrt reflectances', coords={'wavelength':lrt_radf.wavelength})
    refl_lrt.attrs['COT'] = lrt_radf.attrs['COT']
    refl_lrt.attrs['r_eff'] = lrt_radf.attrs['r_eff']                                                                                                                                      
    return refl_lrt;
                          

def calc_HySICS_reflectance(datacube,solar_interp,phase,x,y):    
    f = interp1d(datacube.wavelength,datacube.values[x,y,:],fill_value='NaN')
    hy_new = f(solar_interp.wavelength)
    
    mu = np.cosd(38.5) if phase == 'Liquid Water' else np.cosd(56)
    refl_HySICS= [np.pi/mu*hy_new[w]/solar_interp.values[w] for w in np.arange(0,len(solar_interp.wavelength))]
    
    refl_hysics_pixel = xr.DataArray(refl_HySICS,dims=['wavelength'], name='HySICS reflectance', coords={'wavelength':solar_interp.wavelength})
    return refl_hysics_pixel;


def find_retrieval_wavelengths(num_wl): 
    if num_wl%5 !=0: print('Number of wavelengths in algorithm must be a factor of 5.')   
    
    if (num_wl ==5):
        retrieval_wl_list=[[750], [1000], [1200], [1660], [2200]]
    else:
        wl_lims = [[600,750], [1000,1080], [1240,1320], [1600,1750], [2100,2200]]
        num = num_wl/5
        retrieval_wl_list = [np.linspace(wl_lims[w][0],wl_lims[w][1],num) for w in np.arange(0,5)]
    
    return retrieval_wl_list;


def find_retrieval_idx_list(retrieval_wl_list,solar_interp):   
    idx_list = []
    for i in np.arange(0,5):
        cols = []
        for j in np.arange(0,len(retrieval_wl_list[0])):
            a = [np.min(np.abs(retrieval_wl_list[i][j] - x)) for x in solar_interp.wavelength] 
            val, idx = min((val, idx) for (idx, val) in enumerate(a))
            cols.append(int(idx))
        idx_list.append(cols)
    return idx_list;


def gen_refl_lrt_list(files,phase,solar_interp,phase):
    refl_lrt_list = []
    lrt_rads = [read_LRT(f,phase) for f in files if f.endswith('.dat')]
    lrt_radfs = [filter_LRT_wl(lrt_rad)for lrt_rad in lrt_rads]
    refl_lrts = [calc_LRT_reflectance(r,solar_interp,phase) for r in lrt_radfs]
    return (refl_lrts);


def find_sigma_k_j(refl_lrt,refl_hysics,idx_1,k,j):
    reflectance_term = (5-k)**2*(refl_lrt.values[j]-refl_hysics.values[j])**2
    ratio_term = (k-1)**2*(refl_lrt.values[j]/refl_lrt.values[idx_1]- refl_hysics.values[j]/refl_hysics.values[idx_1])**2
    sigma_k_j = reflectance_term + ratio_term
    return sigma_k_j


def find_sigma_k(refl_lrt,refl_hysics,retrieval_idx_list,idx_1,i):
    k = i+1
    sigma_k = [find_sigma_k_j(refl_lrt,refl_hysics,idx_1,k,j) for j in retrieval_idx_list[i]] 
    return sigma_k


def calc_disagreement(refl_hysics,refl_lrt,retrieval_idx_list):  
    idx_1 = retrieval_idx_list[0][0]
    
    sigma_all = [find_sigma_k(refl_lrt,refl_hysics,retrieval_idx_list,idx_1,i) for i in np.arange(0,5)] 
    
    sigma = np.sum(sigma_all)
    chi_sq = sigma/(len(retrieval_idx_list[0])*60) 
    diff = np.sqrt(chi_sq)*100
    
    disagreement = xr.DataArray(diff)
    disagreement.attrs['COT'] = refl_lrt.attrs['COT']
    disagreement.attrs['r_eff'] = refl_lrt.attrs['r_eff']
    return disagreement;


def calc_best_diff(refl_hysics_pixel,refl_lrt_list, retrieval_idx_list):    
    diff_list = [calc_disagreement(refl_hysics_pixel, r, retrieval_idx_list) for r in refl_lrt_list]
        
    best_diff, idx_diff = min((val, idx) for (idx, val) in enumerate(diff_list))
    return(best_diff);