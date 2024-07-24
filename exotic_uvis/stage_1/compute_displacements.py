import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.stats import norm
from scipy import optimize
from photutils.centroids import centroid_com, centroid_2dg, centroid_quadratic
from exotic_uvis.plotting import plot_exposure
from exotic_uvis.plotting import plot_bkg_stars

def track_bkgstars(obs, bkg_stars, window = 15, verbose_plots = 0, check_all = False, output_dir = None):

    """
    
    Function to compute the x & y displacement of a given background star
    
    """

    # intialize and copy images
    stars_pos, abs_pos = [], []
    images = obs.images.data.copy()

    # iterate over all listed background stars
    for i, pos_init in enumerate(bkg_stars):
        
        # initialize position
        pos = []

        # get window limits
        x0, xf = pos_init[0] - window, pos_init[0] + window
        y0, yf = pos_init[1] - window, pos_init[1] + window

        # iterate over all images
        for image in images:

            # define region around background star
            sub_image = image[y0:yf, x0:xf]
            
            # compute centroid
            x1, y1 = centroid_com(sub_image)

            # append location
            pos.append([x0 + x1, y0 + y1])
        
        rel_pos = np.array(pos) - pos[0]
        
        if check_all:
            plot_exposure([images[0]], scatter_data = [x0 + x1, y0 + y1])

        # save background star location as a function of time
        obs["star{}_disp".format(i)] = (("exp_time", "xy"), rel_pos)
        stars_pos.append(rel_pos)
        abs_pos.append(pos)

    stars_pos = np.array(stars_pos)
    mean_pos = np.mean(stars_pos, axis = 0)

    obs["meanstar_disp"] = (("exp_time", "xy"), mean_pos)
    
    # if true, plot the calculated displacements
    if verbose_plots > 0:
        mean_loc = list(np.mean(abs_pos, axis = 1).transpose())
        plot_bkg_stars(image, obs.exp_time.data, mean_loc, mean_pos, stars_pos, output_dir=output_dir)

    
    return pos


