import csv
import numpy as np
import xarray as xr
from scipy.interpolate import interp1d
import numba


def read_verbose(verbosefile,phase):
    data = [x for x in open(verbosefile,'r')]
    a = data.index('Using new intensity correction, with phase functions\n')
    n=4 if phase=='Liquid Water' else 5
    tau = data[a:][67].split('|')[n].split()[0]
    return(tau);


def read_LRT(file,phase):
    fs = file.split('_')
    reff = fs[3][1:] if hysics_dict['phase']=='Liquid Water' else fs[2][1:]
    ff = file[0:-4] if hysics_dict['phase']=='Liquid Water' else file[0:-9]
    verbosefile = verbose_path+ff+'_550nm_verbose.txt'
    tau = read_verbose(verbosefile,phase)
    
    data = [x for x in csv.reader(open(file,'r'),delimiter='\t')]  
    data = data[5:]
    wl=[]           #Wavelength (nm)
    dirL=[]         #Direct Radiance (W/m^2/nm/sr)
    
    loops = int(len(data)/3)
    for n in range(loops):
        if phase =='Liquid Water':           
            nwl = np.float(data[3*n+1][0][1:9])
            ndirL = np.float(data[3*n+3][0][9:24]) 
        else:
            nwl = np.float(data[3*n][0][1:9])
            ndirL = np.float(data[3*n+2][0][9:24])

        wl.append(nwl)
        dirL.append(ndirL)

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


def interpolate_solar(lrt_rad,solar)
    f = interp1d(solar.wavelength,solar.values,fill_value='NaN')
    solar_new = f(lrt_rad.wavelength)
    solar_interp = xr.DataArray(flux_new,dims=['wavelength'], name='interpolated solar flux', coords={'wavelength':lrt_rad.wavelength})
    return solar_interp;


def calc_LRT_reflectance(lrt_rad,solar_interp):
    refl_lrt = [np.pi*rlrt_rad.values[w]/solar_interp.values[w]for w in np.arange(0,len(lrt_rad.wavelength)]
                                                                        
    refl_lrt = xr.DataArray(refl_lrt,dims=['wavelength'], name='lrt reflectances', coords={'wavelength':lrt_rad.wavelength})
    refl_lrt.attrs['COT'] = lrt_rad.attrs['COT']
    refl_lrt.attrs['r_eff'] = lrt_rad.attrs['r_eff']                                                                                                                                      
    return refl_lrt;
                          

def calc_HySICS_reflectance(datacube,solar_interp,x,y):    
    f = interp1d(datacube.wavelength,datacube.values[x,y,:],fill_value='NaN')
    hy_new = f(solar_interp.wavelength)
    
    refl_HySICS= [np.pi*hy_new[w]/solar_interp.values[w] for w in np.arange(0,len(solar_interp.wavelength))]
    
    refl_hysics_pixel = xr.DataArray(refl_HySICS,dims=['wavelength'], name='HySICS reflectance', coords={'wavelength':solar_interp.wavelength})
    return refl_hysics_pixel;


def find_retrieval_wavelengths(num_wl): 
    if num_wl%5 !=0: print('Number of wavelengths in algorithm must be a factor of 5.')
    wl_lims = [[600,750], [980,1050], [1220,1320], [1500,1750], [2100,2200]]
    retrieval_wl_list = []
    num = num_wl/5
    for w in np.arange(0,5):
        wl = np.linspace(wl_lims[w][0],wl_lims[w][1],num)
        retrieval_wl_list.append(wl)
    if (num_wl ==5):
        retrieval_wl_list=[[750], [1000], [1200], [1660], [2200]] 
    return retrieval_wl_list;


    def calc_disagreement(refl_hysics, refl_lrt, retrieval_wl_list):  
    idx_1 = np.argmin([np.abs(retrieval_wl_list[0][0] - x) for x in refl_hysics.wavelength])
    
    sigma_all = []
    for i in np.arange(0,5):
        k=i+1
        
        sigma_k = []      
        for j in np.arange(0,len(retrieval_wl_list[i])):
            idx_wl = np.argmin([np.abs(retrieval_wl_list[i][j] - x) for x in refl_hysics.wavelength])
            
            reflectance_term = (5-k)**2*(refl_lrt.values[idx_wl]-refl_hysics.values[idx_wl])**2
            ratio_term = (k-1)**2*(refl_lrt.values[idx_wl]/refl_lrt.values[idx_1]- \
                                refl_hysics.values[idx_wl]/refl_hysics.values[idx_1])**2
            
            sigma_k_j = reflectance_term + ratio_term
            sigma_k.append(sigma_k_j)
        sigma_all.append(sigma_k)
        
    sigma = np.sum(sigma_all)
    chi_sq = sigma/(wl_num*12) 
    diff = np.sqrt(chi_sq)*100
    
    disagreement = xr.DataArray(diff)
    disagreement.attrs['COT'] = refl_lrt.attrs['COT']
    disagreement.attrs['r_eff'] = refl_lrt.attrs['r_eff']
    return disagreement;


def gen_refl_lrts_array(path,phase,solar_interp):
    refl_lrt_list = []
    files = [os.path.join(path,f) for f in os.listdir(path)]
    
    lrt_rads = [read_LRT(f,phase) for f in files]
    refl_lrts = [calc_LRT_reflectance(r,solar_interp) for r in lrt_rads]
    refl_lrt_array = da.stack(refl_lrts)
    return (refl_lrt_arrray);

def calc_disagreement_all_LRT(refl_hysics_pixel,refl_lrt_array, retrieval_wl_list):    
    reff_list = [] ; COT_list = [] ; diff_list = []

    for i in np.arange(0,len(refl_lrt_list)):
        diff_dict = calc_disagreement(refl_hysics_pixel, refl_lrt_array[i], retrieval_wl_list)

        reff_list.append(disagreement.attr['r_eff'])
        COT_list.append(disagreement.attr['COT'])   
        diff_list.append(disagreement.values)
        
    val_diff, idx_diff = min((val, idx) for (idx, val) in enumerate(diff_list))

    reff_float = [float(i) for sublist in reff_list for i in sublist]
    COT_float = [float(i) for i in COT_list]
    reff_best = float(reff_float[idx_diff])
    COT_best = float(COT_float[idx_diff]) 