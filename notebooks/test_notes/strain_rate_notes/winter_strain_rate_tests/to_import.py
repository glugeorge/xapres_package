# to_import
# imports everything such that all the test notebooks don't need to list out all the relevant libraries

import sys
import pandas as pd
sys.path.append("../../../../xapres_package/")
import ApRESDefs
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import importlib
import gcsfs
import math
from datetime import datetime

# To run
def reload_winter(site):
    filename = f'gs://ldeo-glaciology/apres/greenland/2022/single_zarrs_noencode/{site}_winter22_23'
    ds = xr.open_dataset(filename,
        engine='zarr', 
        consolidated=True, 
        chunks={}) 
    return ds

def convert_to_seconds(s):
    seconds_per_unit = {"S": 1, "M": 60, "H": 3600, "D": 86400, "W": 604800}
    return int(s[:-1]) * seconds_per_unit[s[-1]]

def sum_error(errs):
    return np.sqrt(np.sum(errs**2))

def error_prop(data):
    return xr.apply_ufunc(sum_error, data, input_core_dims=[["time"]], vectorize = True)

def generate_plots(station,stacking,vd_window,vd_sum,vd_coarse):
    xa = ApRESDefs.xapres(loglevel='debug')
    ds_winter = reload_winter(station)
    current_time = datetime.now().strftime("%H:%M:%S")
    print(f'{current_time}: Data loaded')

    # Determine profiles
    profiles = ds_winter.isel(attenuator_setting_pair=1,time=range(1,len(ds_winter.time))).profile_stacked.resample(time=stacking).mean(dim='time').compute()

    current_time = datetime.now().strftime("%H:%M:%S")
    print(f'{current_time}: Profiles calculated')

    # Calculate displacements
    b1 = profiles.isel(time=range(0,len(profiles.time)-1)).where(profiles.profile_range >= 10,drop=True).compute()
    b2 = profiles.isel(time=range(1,len(profiles.time))).where(profiles.profile_range >= 10,drop=True).compute()
    ds, co, phi = xa.generate_range_diff(b1,b2,vd_window,vd_window,None,0,0.95)

    current_time = datetime.now().strftime("%H:%M:%S")
    print(f'{current_time}: Displacements calculated')

    disp_stack = ds.range_diff.resample(time=vd_sum).sum(dim='time')/convert_to_seconds(vd_sum)*31536000
    strain_polyfit = disp_stack.where(disp_stack.profile_range <= 500,drop=True).polyfit('profile_range',1)
    err_stack = ds.err.resample(time=vd_sum).apply(error_prop)/convert_to_seconds(vd_sum)*31536000
    current_time = datetime.now().strftime("%H:%M:%S")
    print(f'{current_time}: Strain rate calculated')

    # Plotting
    fig0, axs = plt.subplots(nrows=3,figsize=(15,15),sharex=True,layout='constrained')
    fig0.suptitle(f'Station {station}, Atten: {profiles.attenuator.values}dB, Gain: {profiles.AFGain.values}dB, stacked over {stacking}, strain rates over {vd_sum}',fontsize=16)

    profile_to_plot = xa.dB(profiles)
    profile_to_plot.plot(ax=axs[0],x='time',cbar_kwargs={'location':'bottom','fraction':0.05,'pad':0.01,'label':'Reflector amplitude [dB]'})

    axs[0].invert_yaxis()
    axs[0].set_title('Reflector Amplitudes')
    axs[0].set_xlabel('')

    for i in np.arange(0,len(ds.profile_range),len(ds.profile_range)//10):
        (100*ds.range_diff.cumsum(dim='time').isel(profile_range = i) + ds.profile_range.isel(profile_range = i)).plot(ax=axs[1],color='k')
    axs[1].invert_yaxis()
    axs[1].set_title('Reflector Displacements')
    axs[1].set_xlabel('')
    axs[1].set_ylabel('depth [m], displacement [cm]')


    strain_polyfit.polyfit_coefficients.sel(degree=1).plot(ax=axs[2],color='k')
    axs[2].set_title('Strain rate over time')
    axs[2].set_xlabel('')
    axs[2].set_ylabel('strain rate [y$^{-1}$]')

    fig0.supxlabel('date')

    fig1, axs_1 = plt.subplots(ncols=8,figsize=(15,5),sharey=True,layout='constrained')
    fig1.suptitle('Strain fitting')
    # Plotting displacement for individual layers 
    plt_count = 0
    for i in np.arange(math.ceil(0.05*len(disp_stack.time)),len(disp_stack.time),len(disp_stack.time)//8):
        axs_1[plt_count].errorbar(disp_stack.isel(time = i) ,disp_stack.profile_range,yerr=None,xerr=err_stack.isel(time = i),linestyle='',marker='.',color='k')
        fit = strain_polyfit.polyfit_coefficients.sel(degree=1).isel(time=i)*disp_stack.profile_range+strain_polyfit.polyfit_coefficients.sel(degree=0).isel(time=i)

        axs_1[plt_count].plot(fit,disp_stack.profile_range,color='r')
        unit = 'yr$^{-1}$'
        axs_1[plt_count].text(np.percentile(disp_stack.isel(time = i),7),0, f"\u03B5\u0307= {strain_polyfit.polyfit_coefficients.sel(degree=1).isel(time=i).values:.3g}"+unit)
        axs_1[plt_count].set_title(f'{disp_stack.time.isel(time = i).values.astype(str)[:10]}')
        axs_1[plt_count].set_xlim(np.percentile(disp_stack.isel(time = i),[5,95]))
        plt_count += 1
        if plt_count == 8:
            break

    axs_1[0].invert_yaxis()

    axs_1[0].set_ylabel('depth [m]')
    fig1.supxlabel('vertical velocity [m y$^{-1}$]')