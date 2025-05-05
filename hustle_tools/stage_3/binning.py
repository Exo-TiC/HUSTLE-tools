import numpy as np

import xarray as xr


def get_wlc(specs, specs_err, normalize, norm_lim):
    """Computes the normalized white light curve
    using the full 1D spectra.

    Args:
        specs (np.array): spec.data of the obs xarray.
        specs_err (np.array): spec_err.data of the obs xarray.
        normalize (bool, optional): whether to normalize extracted light
        curves. Defaults to True.
        norm_lim (int): how many indices are considered out of transit.

    Returns:
        np.array, np.array: white light curve and squared uncertainties
        normalized by the out-of-transit flux.
    """
    
    # normalization factor
    norm_factor = 1
    if normalize:
        norm_factor = np.median(np.sum(specs[:norm_lim, :], axis = 1))
    
    wlc = np.sum(specs, axis = 1) / norm_factor
    err = np.sqrt(np.sum(specs_err**2, axis = 1)) / norm_factor

    return wlc, err


def get_speclcs(specs, specs_err, waves,
                bin_method, ncol, bins,
                normalize, norm_lim):
    """Generates the spectral light curves from a series of 1D spectra.

    Args:
        specs (np.array): spec.data of the obs xarray.
        specs_err (np.array): spec_err.data of the obs xarray.
        waves (np.array): wave.data of the obs xarray.
        bin_method (str): technique to use for binning. Can be either
        "wavelengths" or "columns".
        ncol (int): if bin_method is "columns", how many detector columns
        go into each bin.
        bins (np.array or float): if bin_method is "wavelengths", either
        the bin edges pre-defined or the spacing between each bin edge,
        to be extrapolated from the lower and upper bounds of waves.
        normalize (bool, optional): whether to normalize extracted light
        curves. Defaults to True.
        norm_lim (int): how many indices are considered out of transit.

    Returns:
        np.array, np.array, np.array, np.array: the binned light curves,
        flux uncertainties, and wavelength centers/bins for each curve.
    """

    lc_binned, lc_error, wave_edges_acc= [], [], []

    waves_mid = (waves[:-1] + waves[1:])/2
    waves_mid =  np.concatenate(([2*waves[0] - waves_mid[0]], waves_mid, [2*waves[-1] - waves_mid[-1]]))

    # check bin technique
    if bin_method == 'wavelengths':
        # build our wavelength bins object, if needed
        if np.isscalar(bins):
            wave_edges = np.arange(waves[0], waves[-1], bins)
        else:
            wave_edges = bins

    elif bin_method == 'columns':
        # build our wavelength bins object by column index
        wave_edges = np.arange(0,len(waves),ncol)
        wave_edges = waves[wave_edges]
        # grab the last wavelength
        wave_edges = np.append(wave_edges,waves[-1])

    # get the actual bin edges and centers    
    wave_cents = (wave_edges[1:] + wave_edges[:-1])/2
    wave_edges_acc.append(waves_mid[0])

    # and bin
    for i in range(len(wave_cents)):

        mask = (waves >= wave_edges[i]) & (waves < wave_edges[i + 1]) # do this more accurately
        if not np.any(mask==True):
            # The entire mask is false, so must pass and delete this wave cent.
            wave_cents = np.delete(wave_cents,i)
            continue

        wave_edges_acc.append(waves_mid[np.where(mask == True)[0][-1] + 1])

        norm_factor = 1
        if normalize:
            norm_factor = np.median(np.sum(specs[:norm_lim, mask], axis = 1))
        lc_binned.append(np.sum(specs[:, mask], axis = 1) / norm_factor)
        lc_error.append(np.sqrt(np.sum(specs_err[:, mask]**2, axis = 1)) / norm_factor)
    
    # array-ify it
    wave_edges_acc = np.array(wave_edges_acc)
    wave_cents_acc = (wave_edges_acc[1:] + wave_edges_acc[:-1])/2

    return np.array(lc_binned), np.array(lc_error), wave_cents, wave_edges


def bin_light_curves(obs, order, bin_method,
                     bins = 100, ncol = 10,
                     rem_exp = None, normalize = True, norm_lim = 10):
    """Generates all white light and spectroscopic light curves according
    to the selected binning techniques and data trimming requests

    Args:
        obs (xarray): contains spec.data, spec_err.data, exp_time.data,
        and wave.data which are all needed for binning.
        order (str): used to ensure that the waves are sorted lowest to
        highest, which they are in negative orders but aren't in positive.
        bin_method (str): technique to use for binning. Can be either
        "wavelengths" or "columns".
        ncol (int): if bin_method is "columns", how many detector columns
        go into each bin.
        bins (np.array or float): if bin_method is "wavelengths", either
        the bin edges pre-defined or the spacing between each bin edge,
        to be extrapolated from the lower and upper bounds of waves.
        rem_exp (np.array, optional): list of exposure times to remove, if
        any frames need to be kicked. Defaults to None.
        normalize (bool, optional): whether to normalize extracted light
        curves. Defaults to True.
        norm_lim (int, optional): how many indices are considered out of transit
        or out of eclipse. Defaults to 10.


    Returns:
        xarray: the broadband and spectroscopic light curves for this order.
    """
    # Remove exposures:
    if rem_exp:
        obs = obs.drop_isel(exp_time = rem_exp)
    
    # load spectra data
    specs = obs.spec.data
    specs_err = obs.spec_err.data
    
    # load time
    exp_times = obs.exp_time.data

    # get wavelengths
    waves = obs.wave.data

    # reverse wavelengths if positive order
    if order[0] == '+':
        sort_vect = np.argsort(waves)
        waves = waves[sort_vect]
        specs = specs[:, sort_vect]
        specs_err = specs_err[:, sort_vect]

    # compute white light curve
    wlc, wlc_err = get_wlc(specs, specs_err, normalize, norm_lim)

    # compute spectroscopic light curves
    spec_lcs, spec_lcs_errs, wave_cents, wave_edges = get_speclcs(specs, specs_err, waves,
                                                                  bin_method, ncol, bins,
                                                                  normalize, norm_lim)
    #print(np.shape(exp_times), np.shape(specs), np.shape(spec_lcs), np.shape(spec_lcs_errs), np.shape(wave_cents), np.shape(wave_edges))

    # create xarray with data:
    light_curves = xr.Dataset(
            data_vars=dict(
                wlc = (['exp_time'], wlc),
                wlc_err = (['exp_time'], wlc_err),
                spec_lc = (['wave_cents', 'exp_time'], spec_lcs),
                spec_err = (['wave_cents', 'exp_time'], spec_lcs_errs),
                ),
            coords=dict(
                exp_time = exp_times,
                wave_cents = wave_cents,
                wave_edges = wave_edges,
        ),
    ) 

    return light_curves

