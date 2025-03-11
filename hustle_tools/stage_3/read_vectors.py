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


def correct_vects(exp_times, vect):


    off1 = vect[26] - vect[25]
    off2 = 0
    off3 = 0
    
    #coeffs = np.polyfit(exp_times[-8:-1], vect[-8:-1], deg = 6)
    #extrapol = np.polyval(coeffs, exp_times[-1])
    
    vect_corr = vect.copy()
    vect_corr[-1] = vect[-2] + off1

    plt.figure(figsize = (10, 7))
    plt.scatter(exp_times, vect)
    #plt.plot(exp_times[-8:], np.polyval(coeffs, exp_times[-8:]))
    plt.scatter(exp_times[-1], vect_corr[-1])
    plt.show()


    return [vect_corr]



def draw_ephemeris(ephemeris, period):

    """
    
    Function to draw the ephemerides and find the correct mid transit time
    
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

    """
    
    Function to normalize a set of arrays
    
    """

    medians = np.median(vects, axis = 1)[:, np.newaxis]
    stds = np.std(vects, axis = 1)[:, np.newaxis]
    state_vects = (vects - medians)/stds

    return state_vects

 

def get_jitter_data(data_dir, res):

    """
    
    Function to get the jitter data
    
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



def get_state_vects(res, data_dir = None, method = 'jit_dec', include = 'all', plot = False):

    """
    
    Function to prepare the state vectors
    
    """
    
    # calculate phase vector
    HST_hperiod = 2853.139318 #what is the exact hst orbit time? (95.1, 95.7, 96?)
    phase_offset = 0.
    HST_orbit = HST_hperiod*2 / (24 * 3600) # HST orbit in days
    hst_phase = ((res.exp_time.data + phase_offset) % HST_orbit) / HST_orbit

    # import image parameters
    #widths = np.mean(res.fit_widths.data * res.spec.data, axis = 1)
    widths = np.mean(res.fit_widths.data, axis = 1)
    bkg_star_x = res.meanstar_disp.data[:, 0]
    bkg_star_y = res.meanstar_disp.data[:, 1]
    spec_disp = res.spec_disp.data #is this correctly flipped?
    prof_disp = np.median(res.prof_disp.data, axis = 1)

    # get jitter vectors
    jitter_vects, jitter_names = get_jitter_data(data_dir, res)

    # basic vectors
    emp_vects = np.array([res.exp_time.data,
                          hst_phase,
                          spec_disp])
    vect_names = ['Time', 'HST phase', 'Wave shift']
    

    # if true, include jitter decorrelation vevtors
    if method == 'jit_dec':

        if include == 'all':
            extra_vects = np.array([
                            widths,
                            bkg_star_x,
                            bkg_star_y])
            extra_vects = np.concatenate((extra_vects, jitter_vects))

        else:
            extra_vects = np.array([
                            hst_phase**2,
                            hst_phase**3,
                            hst_phase**4,
                            widths,
                            #bkg_star_x,
                            #bkg_star_y,
                            prof_disp,
                            ])
            vect_names.extend(['HST phase**2', 'HST phase**3', 'HST phase**4', 'Prof widths', 'Prof shifts'])

            for var in include:
                ind = np.where(jitter_names == var)
                jit_vect = jitter_vects[ind]
                #jit_vect = correct_vects(res.exp_time.data, jit_vect.flatten())
                extra_vects = np.concatenate((extra_vects, jit_vect))
                vect_names.append(var)
        


    elif method == 'sys_marg':
        extra_vects = np.array([
                        hst_phase**2,
                        hst_phase**3,
                        hst_phase**4,
                        widths, 
                        bkg_star_x,
                        bkg_star_y, 
                        prof_disp])  


    all_vects = np.concatenate((emp_vects, extra_vects))
    state_vects = normalize(all_vects)

    
    if plot:
        # plot main detrending vectors:
        plt.figure(figsize = (15, 12))
        for i, state_vect in enumerate(state_vects):
            #state_vects[8 + i, -1] = np.median(state_vect)
            if i < 16:
                ax = plt.subplot(4, 4, i + 1)
                ax.scatter(res.exp_time.data, state_vect)
        plt.show()

    return state_vects, vect_names