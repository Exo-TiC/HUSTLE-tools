import os
from tqdm import tqdm

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from matplotlib.pyplot import rc

rc('font', family='serif')
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('axes', labelsize=14)
plt.rc('legend',**{'fontsize':11})


from hustle_tools.plotting import plot_exposure, plot_0th_histogram, plot_0th_profile


def Gauss1D(x, A, x0, sigma):
    """Creates a 1D Gaussian on the given x range.

    Args:
        x (np.array): independent variable in the Gaussian.
        H (float): vertical offset.
        A (float): amplitude of the Gaussian.
        x0 (float): center of the Gaussian.
        sigma (float): width of the Gaussian.

    Returns:
        np.array: Gaussian profile on domain x.
    """
    return A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def exponential(x, x0, y0, a, b):
    """Creates a 1D exponential on a given x range.

    Args:
        x (np.array): independent variable in the Gaussian.
        x0 (float): normalisation of independent variable.
        y0 (float): constant offset of the exp.
        a (float): amplitude of the exponential.
        b (float): power of the exponential.

    Returns:
        p.array: 1D exponential on domain x.
    """
    return y0 + a * np.exp(b * (x - x0)) 


def remove_zeroth_order(obs, mode = 'radial_profile', zero_pos = [1158, 300],
                        rmin = 100, rmax = 300, rwidth = 3, fit_profile = False, 
                        verbose = 0, show_plots = 0, save_plots = 0, output_dir = None):
    """Method to model and remove the 0th order wings, which can contaminate
    the bluest end of the +1 and -1 orders.

    Args:
        obs (xarray): images to remove the 0th order from.
        mode (str, optional): type of profile to model the 0th with.
        Defaults to 'radial_profile'.
        zero_pos (list, optional): _description_. Defaults to [1158, 300].
        rmin (int, optional): _description_. Defaults to 100.
        rmax (int, optional): _description_. Defaults to 300.
        rwidth (int, optional): _description_. Defaults to 3.
        fit_profile (bool, optional): _description_. Defaults to False.
        verbose (int, optional): How detailed you want the printed statements
        to be. Defaults to 0.
        show_plots (int, optional): How many plots you want to show. Defaults to 0.
        save_plots (int, optional): How many plots you want to save. Defaults to 0.
        output_dir (str, optional): Where to save the plots to, if save_plots
        is greater than 0. Defaults to None.

    Returns:
        array-like: obs corrected and a model of the 0th order in each frame.
    """
    
    if verbose>0:
        print("Removing 0th order contamination...")

    # copy images array
    images = obs.images.data.copy()

    if mode == 'radial_profile':

        # initialize zeroth background image
        zero_bkg = np.zeros_like(images)

        # initialize x and y arrays
        xvals = np.arange(np.shape(images)[2])
        yvals = np.arange(np.shape(images)[1])

        # create polar coordinates grids
        xx, yy = np.meshgrid(xvals, yvals)
        rr = np.sqrt((xx - zero_pos[0])**2 + (yy - zero_pos[1])**2)
        tt = np.arctan2((yy - zero_pos[1]), (xx - zero_pos[0]))

        # define radial vectors
        profiles = []
        r_vect = np.arange(rmin, rmax, rwidth)
        r_cents = (r_vect[:-1] + r_vect[1:])/2

        # define region for estimating zeroth profile
        ang = np.pi/6
        mask_angle = ((tt > ang) & (tt < np.pi - ang)) | ((tt > ang - np.pi) & (tt < -ang))

        if verbose == 2:
            print("Radial profile mask_angle sum:",np.sum(mask_angle))

        # plot the area used to build the histogram
        if (save_plots>0 or show_plots>0):
            hist_area = images[0].copy()
            hist_area[mask_angle & (rr < rmax) & (rr > rmin)] = 'nan'

            plot_exposure([hist_area], title = 'Area to build the histogram', 
                        show_plot=(show_plots>0), save_plot=(save_plots>0),
                        output_dir=output_dir, filename = ['0th_order_histogram-area_max'])
        

        # iterate over all images
        for j, image in enumerate(tqdm(images,
                                       desc = 'Calculating 0th order background... Process:',
                                       disable=(verbose==0))):
            # initialize profile
            profile = []
        
            #if j >= 0:
            # iterate over all radial positions
            for i in range(len(r_vect) - 1):
                
                # define mask of area to compute radial profile
                mask_radius = (rr >= r_vect[i]) & (rr < r_vect[i + 1]) 
                mask = mask_radius & mask_angle & (image != 0) #& (image < 60)

                # create histogram
                hist, bin_edges = np.histogram(image[mask], bins = np.linspace(-10, 60, 100))
                bin_cents = (bin_edges[:-1] + bin_edges[1:])/2

                # do gaussian fit to histogram
                try:
                    parameters, covariance = optimize.curve_fit(Gauss1D, bin_cents, 
                                                                hist, 
                                                                p0 = [np.amax(hist), bin_cents[np.argmax(hist)], 10])
                    rad_value = parameters[1]
                    gaussian_curve = Gauss1D(bin_cents, parameters[0], parameters[1], parameters[2])
                except:
                    gaussian_curve=None
                    rad_value = np.median(image[mask])
                
                # plot one case half way through
                if i == int(len(r_vect)/2):
                    hist_area = image.copy()
                    hist_area[mask_angle & (rr < r_vect[i + 1]) & (rr > r_vect[i])] = 'nan'

                    if (save_plots==2 or show_plots==2):
                        plot_exposure([hist_area], title = 'Area to build the histogram', 
                                    show_plot=(show_plots==2), save_plot=(save_plots==2),
                                    output_dir=output_dir, filename = [f'0th_order_histogram-area_frame{j}'])

                        plot_0th_histogram(image[mask], bin_cents, gaussian_curve, rad_value, 
                                            show_plot = (show_plots > 1), save_plot = (save_plots > 1),
                                            filename = f"0th_order_histogram_frame{j}.png", output_dir = output_dir)
                        
                profile.append(rad_value)

                if fit_profile == False:
                    image[mask_radius] -= rad_value 
                    zero_bkg[j, mask_radius] = rad_value

            # fit radial profile with exponential for smoothness
            if fit_profile:
                
                # perform curve fit on extracted radial profile
                profile = np.array(profile)
                mask_radius = (rr >= rmin) & (rr < rmax) 
                mask = r_cents > 0
                popt, pcov = optimize.curve_fit(exponential, r_cents[mask],
                                                profile[mask], p0 = [100., 7., 1., -0.1])
                
                # remove profile from exposures
                fitted_profile = exponential(r_cents, popt[0], popt[1], popt[2], popt[3])
                image[mask_radius] -= exponential(rr[mask_radius], popt[0], popt[1], popt[2], popt[3])
                zero_bkg[j, mask_radius] = exponential(rr[mask_radius], popt[0], popt[1], popt[2], popt[3])

                # plot extracted profile
                if (save_plots==2 or show_plots==2):
                    plot_0th_profile(r_cents, profile, fitted_profile=None,
                                        show_plot = (show_plots > 1), save_plot = (save_plots > 1),
                                        filename =f"0th_order_radial-vals_frame{j}.png", output_dir = output_dir)

            profiles.append(profile)
            
            # do basic plots
            if (j == 0) and ((show_plots>0) or (save_plots>0)):
                plot_exposure([obs.images.data[j, :], image], title = '0th order removal example', 
                                show_plot=(show_plots>0), save_plot=(save_plots>0),
                                output_dir=output_dir, filename = ['0th_order_corrected_frame0'])
                
                plot_exposure([zero_bkg[j]], title = '0th order model', #min=np.min(fitted_profile), max=np.max(fitted_profile),
                                show_plot=(show_plots>0), save_plot=(save_plots>0),
                                output_dir=output_dir, filename = ['0th_order_model_frame0'])
                
            # do additional plots
            if ((show_plots==2) or (save_plots==2)):
                plot_exposure([obs.images.data[j, :], image], title = '0th order removal, frame {}'.format(j), 
                                show_plot=(show_plots==2), save_plot=(save_plots==2),
                                output_dir=output_dir, filename = [f'0th_order_corrected_frame{j}'])
                
                plot_exposure([zero_bkg[j]], title = '0th order model', #min=np.min(fitted_profile), max=np.max(fitted_profile),
                                show_plot=(show_plots>0), save_plot=(save_plots>0),
                                output_dir=output_dir, filename = [f'0th_order_model_frame{j}'])
        
        # transform to numpy array and save data            
        profiles = np.array(profiles)
        obs.images.data = images

    return zero_bkg
