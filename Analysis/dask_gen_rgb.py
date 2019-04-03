import xarray as xr

datacube_dvc = xr.open_dataarray('datacube_dvc.nc',chunks=10)

def two_percent_linear_stretch(x,idx):
    c = x.isel(wavelength=idx)/np.amax(x.isel(wavelength=idx))
    p2, p98 = np.percentile(c, (2, 98))
    c_scaled = exposure.rescale_intensity(c, in_range=(p2, p98))
    return c_scaled;


def find_rgb_idx(wl_list):
    r_idx = np.argmin([np.abs(670 - x) for x in wl_list])
    g_idx = np.argmin([np.abs(555 - x) for x in wl_list])
    b_idx = np.argmin([np.abs(443 - x) for x in wl_list])
    rgb_idx_list = [r_idx,g_idx,b_idx]
    return (rgb_idx_list); 


def apply_stretch(datacube):
    rgb_idx_list = find_rgb_idx(datacube.coords['wavelength'])
    
    red = two_percent_linear_stretch(datacube,rgb_idx_list[0])
    green = two_percent_linear_stretch(datacube,rgb_idx_list[1])
    blue = two_percent_linear_stretch(datacube,rgb_idx_list[2])

    rgb = np.dstack((red,green,blue))
    rgb = xr.DataArray(rgb,dims=('x','y','rgb'))
    return rgb;

rgb = apply_stretch(datacube_dvc)
rgb.to_netcdf('rgb_dvc.nc',format='NETCDF4')