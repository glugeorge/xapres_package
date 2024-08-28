import sys
import pandas as pd
sys.path.append("../../../xapres_package/")
import ApRESDefs
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import importlib
import gcsfs
import math
from datetime import datetime

def reload(site):
    filename = f'gs://ldeo-glaciology/apres/greenland/2022/single_zarrs_noencode/{site}'
    ds = xr.open_dataset(filename,
        engine='zarr', 
        chunks={}) 
    return ds
    
def reload_winter(site):
    filename = f'gs://ldeo-glaciology/apres/greenland/2022/single_zarrs_noencode/{site}_winter22_23'
    ds = xr.open_dataset(filename,
        engine='zarr', 
        consolidated=True, 
        chunks={}) 
    return ds

def reload_summer(site):
    filename = f'gs://ldeo-glaciology/apres/greenland/2022/single_zarrs_noencode/{site}_summer_23'
    ds = xr.open_dataset(filename,
        engine='zarr', 
        consolidated=True, 
        chunks={}) 
    return ds

def plot_hist(ds,ylim=False):
    atten_count = len(ds.attenuator_setting_pair)
    fig, axs = plt.subplots(ncols=atten_count,figsize=(5*atten_count,5))
    
    for i in range(atten_count):
        ds.chirp.isel(attenuator_setting_pair=i).mean(dim='chirp_num').plot.hist(ax=axs[i],bins=100)
        axs[i].set_title(f'G = {int(ds.AFGain[i].values)} dB, A = {int(ds.attenuator[i].values)} dB')
        axs[i].set_ylabel('count')
        axs[i].set_xlabel('chirp voltage [V]')
        if ylim:
            old_lim = axs[i].get_ylim()
            axs[i].set_ylim([old_lim[0],old_lim[1]*0.2])

def plot_amplitude_trends(ds):
    atten_count = len(ds.attenuator_setting_pair)
    fig, axs = plt.subplots(nrows=3,ncols=atten_count,figsize=(5*atten_count,15))
    for i in range(atten_count):
        # Over chirp 
        abs(ds.isel(attenuator_setting_pair = i).chirp).max(dim=['chirp_num','time'])[1:-1].plot(ax=axs[0] [i],linestyle='None',marker='.',label='max')
        abs(ds.isel(attenuator_setting_pair = i).chirp).mean(dim=['chirp_num','time'])[1:-1].plot(ax=axs[0][i],linestyle='None',marker='.',label='mean')
        axs[0][i].set_title(f'G = {int(ds.AFGain[i].values)} dB, A = {int(ds.attenuator[i].values)} dB')
        axs[0][i].set_ylabel('')
        axs[0][i].set_xlabel('chirp time')
        axs[0][i].legend()

        # Over chirp number 
        abs(ds.isel(attenuator_setting_pair = i).chirp).max(dim=['chirp_time','time']).plot(ax=axs[1][i],linestyle='None',marker='.',label='max')
        abs(ds.isel(attenuator_setting_pair = i).chirp).mean(dim=['chirp_time','time']).plot(ax=axs[1][i],linestyle='None',marker='.',label='mean')
        axs[1][i].set_ylabel('')
        axs[1][i].set_xlabel('chirp number')
        axs[1][i].set_title('')

        # over burst time
        abs(ds.isel(attenuator_setting_pair = i).chirp).max(dim=['chirp_time','chirp_num']).plot(ax=axs[2][i],linestyle='None',marker='.',label='max')
        abs(ds.isel(attenuator_setting_pair = i).chirp).mean(dim=['chirp_time','chirp_num']).plot(ax=axs[2][i],linestyle='None',marker='.',label='mean')
        axs[2][i].set_ylabel('')
        axs[2][i].set_xlabel('date')
        axs[2][i].set_title('')
    
    axs[0][0].set_ylabel('chirp magnitude over chirp time [V]')
    axs[1][0].set_ylabel('chirp magnitude over chirp number [V]')
    axs[2][0].set_ylabel('burst magnitude over time [V]')

def plot_bad_chirp_count(ds):
    # input ds should be the form of ds.isel(attenuator_setting_pair=i).chirp.where(condition)
    fig, axs = plt.subplots(ncols=3,figsize=(15,5),layout='tight')
    ds.count(dim=['chirp_num','time']).plot(ax=axs[0],linestyle='None',marker='.')
    axs[0].set_ylabel('count')
    axs[0].set_title('Total instances of clipping in chirp section')
    
    (ds.max(dim='chirp_time').count(dim='time')/(ds.max(dim='chirp_time').count())).plot(ax=axs[1],linestyle='None',marker='.')
    axs[1].set_ylabel('fraction of total clipped chirps')
    axs[1].set_title('Clipped chirp distribution by chirp number')
    
    ds.max(dim='chirp_time').count(dim='chirp_num').plot(ax=axs[2],linestyle='None',marker='.')
    axs[2].set_ylabel('count')
    axs[2].set_title('Number of chirps with clipping in each burst')

def custom_profile(chirps,clip_threshold=1.2,min_chirps = 20,start=0,stop=39999,pad=2):
    # This stop of 39999 by default cuts off the last two samples
    # This function requires starting on an even number
    times = chirps.chirp_time.values.astype('float64')/1e9
    regular_freq_range = 2e8+2e8*times
    start_freq = regular_freq_range[start]
    stop_freq = regular_freq_range[stop-1]
    B = stop_freq - start_freq
    CentreFreq = (stop_freq + start_freq)/2
    K = 2e8 # determined from step-up freq (5000 Hz) and step-up time (2.5e-5)
    c0 = 3e8 # speed of light in vaccuum
    ER_ICE = 3.18
    T0 = times[start]
    T1 = times[stop-1]

    # Drop bad bursts
    chirps = chirps.isel(chirp_time=range(start,stop))
    bad_chirps =  chirps.where(abs(chirps) > clip_threshold)
    good_bursts = bad_chirps.max(dim='chirp_time').count(dim='chirp_num') <= 20-min_chirps
    chirps = chirps.where(good_bursts)
    chirps = chirps.where(abs(chirps).max(dim='chirp_time')<clip_threshold)

    # De-mean and detrend
    chirps = chirps - chirps.mean(dim='chirp_time')
    p = chirps.polyfit('chirp_time',1)
    fit = xr.polyval(chirps.chirp_time, p.polyfit_coefficients)
    chirps = chirps - fit
    chirp_stack = chirps.mean(dim='chirp_num',skipna=True)
    
    window = np.blackman(len(chirp_stack.chirp_time))
    win_chirps = chirp_stack*window
    Nt = len(chirp_stack.chirp_time)
    Nt = math.floor(Nt/2) * 2
    Nfft = math.floor(Nt*pad)
    bin2m = c0/(2.*(T1-T0)*pad*math.sqrt(ER_ICE)*K)
    profile_range = np.asarray([i for i in range(Nfft)]) * bin2m      
    profile_range = profile_range[0:math.floor(Nfft/2)-1]
    padchirp = np.zeros((len(chirp_stack.time),Nfft))
    padchirp[:,0:math.floor(Nt/2)] = win_chirps.data[:,math.floor(Nt/2):-1]
    padchirp[:,-math.floor(Nt/2):] = win_chirps.data[:,0:math.floor(Nt/2)]
    p = np.fft.fft(padchirp,axis=1)/Nfft * math.sqrt(2*pad)
    profile = p[:,0:math.floor(Nfft/2)-1]
    m = np.asarray([i for i in range(profile.shape[1])])/pad
    phiref = 2*math.pi*CentreFreq*m/B - m * m * 2*math.pi * K/2/B**2
    profile_ref = profile * np.exp(phiref[np.newaxis,:]*(-1j))
    profile_range = np.asarray([i for i in range(Nfft)]) * bin2m      
    profile_range = profile_range[0:math.floor(Nfft/2)-1]
    n = np.argmin(profile_range<=1200)
    Range = profile_range[:n]
    Profile = profile_ref[:,:n]
    da = xr.DataArray(Profile,
                  dims=['time','profile_range'],
                  coords={'profile_range': Range,
                         'time': chirps.time.data})
    da.attrs['centre_freq']= CentreFreq
    da.attrs['bandwidth']= B
    da.attrs['start_freq']= start_freq
    da.attrs['stop_fre']= stop_freq
    return da