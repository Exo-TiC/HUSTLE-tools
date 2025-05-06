import numpy as np
from scipy.ndimage import median_filter


def clip_light_curves(light_curves, sigma = 5, verbose = 0):
    """Iterates over the time series of each light curve and clips
    outliers at the specified level.

    Args:
        light_curves (xarray): contains the wlc and spec_lc data vars.
        sigma (int, optional): threshold at which to clip an outlier.
        Defaults to 5.
        verbose (int, optional): how detailed you want the printed statements
        to be. Defaults to 0.

    Returns:
        xarray: light curves with outliers masked by running median.
    """
    # Track all outliers clipped.
    n_clipped_total = 0

    # First, clip outliers from the wlc.
    light_curves.wlc.data, n_clipped = clip_one_curve(light_curves.wlc.data,
                                                      n=int(0.2*len(light_curves.wlc.data)),
                                                      sigma=sigma)    
    n_clipped_total += n_clipped

    if verbose == 2:
        print("{} outliers clipped from broadband light curve.".format(n_clipped))
    
    # Then go through each spectroscopic lightcurve.
    for i in range(light_curves.spec_lc.data.shape[0]):
        light_curves.spec_lc.data[i,:], n_clipped = clip_one_curve(light_curves.spec_lc.data[i,:],
                                                                   n=int(0.2*len(light_curves.spec_lc.data[i,:])),
                                                                   sigma=sigma)
        n_clipped_total += n_clipped
        if verbose == 2:
            print("{} outliers clipped from light curve of wavelength {:.2f} AA.".format(n_clipped,
                                                                                         light_curves.wave_cents[i]))
    
    if verbose > 0:
        print("Sigma clipping removed {} outliers from all time series.".format(n_clipped_total))

    return light_curves


def clip_one_curve(light_curve, n=5, sigma=5):
    """Simple function to replace outliers.

    Args:
        light_curve (np.array): spectral time series.
        n (int, optional): window size. Defaults to 5.
        sigma (int, optional): threshold at which to clip an outlier.
        Defaults to 5.
    
    Returns:
        np.array: time series with outliers masked.
    """
    smoothed = median_filter(light_curve, size=n, mode='wrap')
    n_clipped = np.count_nonzero(np.where(np.abs(smoothed-light_curve)>np.std(light_curve)*sigma,
                                          1,0))
    corrected_curve = np.where(np.abs(smoothed-light_curve)>np.std(light_curve)*sigma,
                               smoothed,light_curve)
    return corrected_curve, n_clipped