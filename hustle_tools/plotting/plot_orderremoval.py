import os
import numpy as np
import matplotlib.pyplot as plt

#define plotting parameters
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('axes', labelsize=14)
plt.rc('legend',**{'fontsize':11})


def plot_0th_histogram(hist_vals, bin_cents, gaussian_fit, fit_value, 
                        show_plot = False, save_plot = False,
                      filename = None, output_dir = None):
    """Plots the histogram corresponding to the values 
    used to calculate the 0th order radial profile

    Args:
        hist_vals (np.array): values used to compute the histogram
        bin_cents (np.array): histogram bin centers
        gaussian_fit (np.array): fitted gaussian curve
        fit_value (float): radial profile fitted value
        show_plot (bool, optional): whether to interrupt execution to
        show the user the plot. Defaults to False.
        save_plot (bool, optional): whether to save this plot to a file.
        Defaults to False.
        filename (str, optional): name to give this file, if saving.
        Defaults to None.
        output_dir (str, optional): where to save the file, if saving.
        Defaults to None.
    """

    plt.figure(figsize = (10, 7))
    plt.hist(hist_vals, bins = bin_cents, color = 'indianred', alpha = 0.7)
    plt.plot(bin_cents, gaussian_fit, color='gray')
    plt.axvline(np.median(hist_vals), color= 'indianred')
    plt.axvline(fit_value, color= 'gray')
    plt.xlabel('Photons')
    plt.ylabel('Counts')

    if save_plot:
        plot_dir = os.path.join(output_dir,'plots')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        plt.savefig(os.path.join(plot_dir, filename),
                    dpi=300,bbox_inches='tight')

    if show_plot:
        plt.show(block=True)
    
    plt.close() # save memory

    return


def plot_0th_profile(r_cents, profile, fitted_profile=None,
                     show_plot = False, save_plot = False,
                      filename = None, output_dir = None):
    """Plots the radial profile corresponding to the 0th order

    Args:
        r_cents (np.array): radial values
        profile (np.array): 0th order radial profile
        fitted_profile (np.array, optional): Fitted exponential 
        to radial values. Defaults to None.
        show_plot (bool, optional): whether to interrupt execution to
        show the user the plot. Defaults to False.
        save_plot (bool, optional): whether to save this plot to a file.
        Defaults to False.
        filename (str, optional): name to give this file, if saving.
        Defaults to None.
        output_dir (str, optional): where to save the file, if saving.
        Defaults to None.
    """

    plt.figure(figsize=(10, 7))
    plt.plot(r_cents, profile, '-o', color='indianred')

    if fitted_profile is not None:
        plt.plot(r_cents, fitted_profile, color = 'gray')

    plt.xlabel('Distance from 0th order')
    plt.ylabel('0th order value')

    if save_plot:
        plot_dir = os.path.join(output_dir,'plots')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        plt.savefig(os.path.join(plot_dir, filename),
                    dpi=300,bbox_inches='tight')

    if show_plot:
        plt.show(block=True)
    
    plt.close() # save memory

    return 