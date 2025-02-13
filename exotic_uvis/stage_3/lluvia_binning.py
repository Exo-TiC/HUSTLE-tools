import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.stats import norm
from scipy import optimize, signal, ndimage
from scipy.interpolate import interp1d
from exotic_ld import StellarLimbDarkening
from itertools import combinations

import xarray as xr
from scipy.optimize import least_squares, curve_fit
import batman
import time
import os



def get_wlc(specs, specs_err, norm_lim):

    '''
    
    Function to compute the white light curve

    '''
    
    # normalization factor
    norm_factor = np.median(np.sum(specs[:norm_lim, :], axis = 1))

    wlc = np.sum(specs, axis = 1) / norm_factor
    err = np.sqrt(np.sum(specs_err**2, axis = 1)) / norm_factor

    return wlc, err



def get_speclcs(specs, specs_err, waves, bins, norm_lim = 8):

    '''
    
    Function to generate the spectral light curves
    
    '''

    lc_binned, lc_error, wave_edges_acc= [], [], []

    waves_mid = (waves[:-1] + waves[1:])/2
    waves_mid =  np.concatenate(([2*waves[0] - waves_mid[0]], waves_mid, [2*waves[-1] - waves_mid[-1]]))

    if np.isscalar(bins):
        wave_edges = np.arange(waves[0], waves[-1], bins)
    else:
        wave_edges = bins
        
    wave_cents = (wave_edges[1:] + wave_edges[:-1])/2
    wave_edges_acc.append(waves_mid[0])

    for i in range(len(wave_cents)):

        mask = (waves >= wave_edges[i]) & (waves < wave_edges[i + 1]) # do this more accurately

        wave_edges_acc.append(waves_mid[np.where(mask == True)[0][-1] + 1])
    
        norm_factor = np.median(np.sum(specs[:norm_lim, mask], axis = 1))
        lc_binned.append(np.sum(specs[:, mask], axis = 1) / norm_factor)
        lc_error.append(np.sqrt(np.sum(specs_err[:, mask]**2, axis = 1)) / norm_factor)
    
    wave_edges_acc = np.array(wave_edges_acc)
    wave_cents_acc = (wave_edges_acc[1:] + wave_edges_acc[:-1])/2

    return np.array(lc_binned), np.array(lc_error), wave_cents, wave_edges



def prepare_lcs(obs, system_params, use_norm = False, rem_exp = None, norm_lim = 10, bins = 100):

    """
    
    Function to prepare the light curves for detrending

    """

    # Remove exposures:
    if rem_exp:
        obs = obs.drop_isel(exp_time = rem_exp)
    
    # Choose spectra to use
    if use_norm:
        specs = obs.norm_spec.data
        specs_err = obs.norm_err.data
    
    else:
        specs = obs.spec.data
        specs_err = obs.spec_err.data
    
   
    # normalize time
    epoch = (np.mean(obs.exp_time.data) - system_params['eph']) / system_params['per']
    ref_time = system_params['per'] * round(epoch) + system_params['eph'] 
    exp_times = obs.exp_time.data - ref_time

    print(ref_time)

    # compute white light curve
    waves = obs.wave.data
    wlc, wlc_err = get_wlc(specs, specs_err, norm_lim)

    '''
    for i, lc in enumerate(np.transpose(specs)):
        print(waves[i], i)
        plt.figure(figsize = (10, 7))
        plt.scatter(obs.exp_time.data,lc)
        plt.show()

    '''
    
    # compute spectroscopic light curves
    spec_lcs, spec_errs, wave_cents, wave_edges = get_speclcs(specs, specs_err, waves, bins)
    
    # create xarray with data:
    light_curves = xr.Dataset(
            data_vars=dict(
                wlc = (['exp_time'], wlc),
                wlc_err = (['exp_time'], wlc_err),
                spec_lc = (['wave_cents', 'exp_time'], spec_lcs),
                spec_err = (['wave_cents', 'exp_time'], spec_errs),
                ),
            coords=dict(
                exp_time = exp_times,
                wave_cents = wave_cents,
                wave_edges = wave_edges,
        ),
        attrs = dict(ref_time=ref_time),
    ) 
    return light_curves

