import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.stats import norm
from scipy import optimize, signal, ndimage
from scipy.interpolate import interp1d
from matplotlib.pyplot import rc
import xarray as xr
from itertools import combinations
import time
import os


def draw_ephemeris(ephemeris, period):
    """Function to draw the ephemerides and find the correct mid transit time

    Args:
        ephemeris (_type_): _description_
        period (_type_): _description_
    
    Returns:
        float: 
    """

    def eph_func(x, t0):
        return t0 + x*period

    times = ephemeris[:, 0]
    errors = ephemeris[:, 1]

    epochs = np.round((times - np.amin(times))/period)

    epochs_hr = np.arange(epochs[-1])

    popt, pcov = optimize.curve_fit(eph_func, epochs, times, sigma = errors)
    print(popt)

    plt.figure()
    #plt.plot(epochs_hr, eph_func(epochs_hr, popt[0]))
    plt.axhline(0)
    plt.errorbar(epochs, (times - eph_func(epochs, popt[0]))*24*60, yerr = errors*24*60, fmt = 'o')
    plt.show()


    return popt[0]


def normalize(vects):
    """Normalizes the state vectors.

    Args:
        vects (np.array): pointing state vectors used for jitter decorrelation.

    Returns:
        np.array: normalized pointing state vectors.
    """
    if len(np.shape(vects)) > 1:
        medians = np.median(vects, axis = 1)[:, np.newaxis]
        stds = np.std(vects, axis = 1)[:, np.newaxis]
        state_vects = (vects - medians)/stds
    else:
        median = np.median(vects)
        std = np.std(vects)
        state_vects = (vects - median)/std


    return state_vects

 
def get_jitter_data(data_dir, res):
    """Function to get the jitter data

    Args:
        data_dir (_type_): _description_
        res (_type_): _description_

    Returns:
        _type_: _description_
    """

    jitter_vals = []
    var_names = []
    first = True

    for filename in np.sort(os.listdir(data_dir)):

        if filename[-9:] == '_jit.fits':
            hdul = fits.open(os.path.join(data_dir, filename))

            #print(hdul.info())
            #print(repr(hdul[1].header))
          
            if hdul[0].header['NEXTEND'] > 1:
         
                for i in range(1, len(hdul)):

                    image = np.array(list(hdul[i].data))
                    mean_image = np.mean(image, axis = 0)
                    jitter_vals.append(mean_image[:22])
                
                if first:
                    for j in range(len(mean_image)):
                        var_names.append(hdul[1].header['TTYPE{}'.format(j+1)])
                    first = False

    jitter_vals = np.array(jitter_vals).transpose()
    jitter_names = np.array(var_names)

    #for i, jit_vec in enumerate(jitter_vals):
    #    plt.figure(figsize = (10, 7))
    #    plt.scatter(range(len(jit_vec)), jit_vec)
    #    plt.title('{}'.format(jitter_names[i]))
    #    plt.show()

    return jitter_vals, jitter_names


def get_state_vectors(res, data_dir = None, method = 'jit_dec', include_jitter = False, plot = False):
    """Function to prepare the state vectors

    Args:
        res (_type_): _description_
        data_dir (_type_, optional): _description_. Defaults to None.
        method (str, optional): _description_. Defaults to 'jit_dec'.
        include_jitter (bool, optional): _description_. Defaults to False.
        plot (bool, optional): _description_. Defaults to False.

    Returns:
        _type_: _description_
    """
    
    # calculate phase vector
    HST_hperiod = 2853.139318 #what is the exact hst orbit time? (95.1, 95.7, 96?)
    phase_offset = 0.
    HST_orbit = HST_hperiod*2 / (24 * 3600) # HST orbit in days
    hst_phase = ((res.exp_time.data + phase_offset) % HST_orbit) / HST_orbit

    # import image parameters
    #widths = np.mean(res.fit_widths.data, axis = 1)
    bkg_star_x = res.meanstar_disp.data[:, 0]
    bkg_star_y = res.meanstar_disp.data[:, 1]
    spec_disp = res.spec_disp.data 
    prof_disp = np.median(res.prof_disp.data, axis = 1)

    # get jitter vectors
    #jitter_vects, jitter_names = get_jitter_data(data_dir, res)

    # if true, include jitter decorrelation vevtors
    #if include_jitter:
    #    for var in include:
    #        ind = np.where(jitter_names == var)
    #        jit_vect = jitter_vects[ind]
            #jit_vect = correct_vects(res.exp_time.data, jit_vect.flatten())
    #        extra_vects = np.concatenate((extra_vects, jit_vect))
    #        vect_names.append(var)

    # create xarray with data:
    state_vectors = xr.Dataset(
            data_vars=dict(
                times = (['exp_time'], normalize(res.exp_time.data)),
                phase = (['exp_time'], normalize(hst_phase)),
                phase2 = (['exp_time'], normalize(hst_phase**2)),
                phase3 = (['exp_time'], normalize(hst_phase**3)),
                phase4 = (['exp_time'], normalize(hst_phase**4)),
                bkg_star_x = (['exp_time'], normalize(bkg_star_x)),
                bkg_star_y = (['exp_time'], normalize(bkg_star_y)),
                spec_disp = (['exp_time'], normalize(spec_disp)),
                prof_disp = (['exp_time'], normalize(prof_disp)),

                ),
            coords=dict(
                exp_time = res.exp_time.data,
        ),
        #attrs = dict(ref_time=ref_time),
    ) 

    if plot:
        # plot main detrending vectors:
        plt.figure(figsize = (15, 12))
        i = 0
        for varname, vect in state_vectors.data_vars.items():
            ax = plt.subplot(3, 3, i + 1)
            ax.plot(res.exp_time.data, vect, 'o', color='indianred', markeredgecolor='black')
            ax.set_title(varname)
            i += 1
        plt.show(block=True)


    return state_vectors
