import numpy as np
import xarray as xr

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
        if not np.any(mask==True):
            # The entire mask is false, so must pass and delete this wave cent.
            wave_cents = np.delete(wave_cents,i)
            continue

        wave_edges_acc.append(waves_mid[np.where(mask == True)[0][-1] + 1])
    
        norm_factor = np.median(np.sum(specs[:norm_lim, mask], axis = 1))
        lc_binned.append(np.sum(specs[:, mask], axis = 1) / norm_factor)
        lc_error.append(np.sqrt(np.sum(specs_err[:, mask]**2, axis = 1)) / norm_factor)
    
    wave_edges_acc = np.array(wave_edges_acc)
    wave_cents_acc = (wave_edges_acc[1:] + wave_edges_acc[:-1])/2

    return np.array(lc_binned), np.array(lc_error), wave_cents, wave_edges



def bin_light_curves(obs, order, rem_exp = None, norm_lim = 10, bins = 100):
    """_summary_

    Args:
        obs (_type_): _description_
        system_params (_type_): _description_
        use_norm (bool, optional): _description_. Defaults to False.
        rem_exp (_type_, optional): _description_. Defaults to None.
        norm_lim (int, optional): _description_. Defaults to 10.
        bins (int, optional): _description_. Defaults to 100.

    Returns:
        _type_: _description_
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
    wlc, wlc_err = get_wlc(specs, specs_err, norm_lim)

    # compute spectroscopic light curves
    spec_lcs, spec_lcs_errs, wave_cents, wave_edges = get_speclcs(specs, specs_err, waves, bins, norm_lim=norm_lim)
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
        #attrs = dict(ref_time=ref_time),
    ) 

    return light_curves