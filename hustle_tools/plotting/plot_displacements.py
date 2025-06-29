import os

import numpy as np
import matplotlib.pyplot as plt

from hustle_tools.plotting import plot_exposure


#define plotting parameters
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('axes', labelsize=14)
plt.rc('legend',**{'fontsize':11})


def plot_bkg_stars(image, exp_time, mean_loc, mean_pos, stars_pos, 
                   show_plot = False, save_plot = False, output_dir = None):
    """Function to plot the displacements of the background stars.

    Args:
        image (array-like): 2D image array showing all of the stars in context.
        exp_times (array-like): BJD exposure times for each frame.
        mean_loc (array-like): mean absolute location of each star.
        mean_pos (array-like): mean relative position of each star.
        stars_pos (array-like): relative position of each star over time.
        show_plot (bool, optional): whether to interrupt execution to
        show the user the plot. Defaults to False.
        save_plot (bool, optional): whether to save the plot to a file.
        Defaults to False.
        output_dir (str, optional): output directory where the plot will be
        saved. Defaults to None.
    """

    # first just plot the example frame and where the stars we tracked are
    plot_exposure([image], scatter_data = mean_loc, title = 'Location of background stars',  
                  show_plot=show_plot, save_plot=save_plot,
                  output_dir=output_dir, filename = ['bkg_stars_location'])
      
    # then we plot the displacements of each star in x and y over time

    # initialize figure and plot x motion of each star
    plt.figure(figsize = (10, 7))
    plt.plot(exp_time, mean_pos[:, 0], '-o')
    plt.plot(exp_time, np.transpose(stars_pos[:, :, 0]), '-o', alpha = 0.5)
    plt.xlabel('Time Of Exposure (MJD)')
    plt.ylabel('X-Pixel Displacement (pixels)')
    plt.title("Background Stars X Drift Over Time")

    if save_plot:
        plot_dir = os.path.join(output_dir, 'plots') 
        filedir = os.path.join(plot_dir, 'bkg_stars_x_displacement.png')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        plt.savefig(filedir, bbox_inches='tight', dpi=300)

    if show_plot:
        plt.show(block=True)

    plt.close() # save memory

    # initialize figure and plot y motion of each star
    plt.figure(figsize = (10, 7))
    plt.plot(exp_time, mean_pos[:, 1], '-o')
    plt.plot(exp_time, np.transpose(stars_pos[:, :, 1]), '-o', alpha = 0.5)
    plt.xlabel('Time Of Exposure (MJD)')
    plt.ylabel('Y-Pixel Displacement (pixels)')
    plt.title("Background Stars Y Drift Over Time")

    if save_plot:
        plot_dir = os.path.join(output_dir, 'plots') 
        filedir = os.path.join(plot_dir, 'bkg_stars_y_displacement.png')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir) 
        plt.savefig(filedir, bbox_inches='tight', dpi=300)

    if show_plot:
        plt.show(block=True)

    plt.close() # save memory
    
    return 


def plot_0th_order(exp_time, xs, ys,
                   show_plot = False, save_plot = False, output_dir = None):
    """Function to plot the location of the 0th order over time.

    Args:
        exp_times (array-like): BJD exposure times for each frame.
        xs (array-like): 0th order x positions.
        ys (array-like): 0th order y positions.
        show_plot (bool, optional): whether to interrupt execution to
        show the user the plot. Defaults to False.
        save_plot (bool, optional): whether to save the plot to a file.
        Defaults to False.
        output_dir (str, optional): output directory where the plot will be
        saved. Defaults to None.
    """

    # initialize figure and plot x motion of 0th order
    plt.figure(figsize = (10, 7))
    plt.plot(exp_time, xs, '-o')
    plt.xlabel('Exposure Times (MJD)')
    plt.ylabel('X-Pixel Displacement (pixels)')
    plt.title("0th Order X Drift Over Time")

    if save_plot:
        plot_dir = os.path.join(output_dir, 'plots') 
        filedir = os.path.join(plot_dir, '0th_order_x_displacement.png')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir) 
        plt.savefig(filedir, bbox_inches='tight', dpi=300)

    if show_plot:
        plt.show(block=True)

    plt.close() # save memory

    # initialize figure and plot y motion of 0th order
    plt.figure(figsize = (10, 7))
    plt.plot(exp_time, ys, '-o')
    plt.xlabel('Exposure Times (MJD)')
    plt.ylabel('Y-Pixel Displacement (pixels)')
    plt.title("0th Order Y Drift Over Time")

    if save_plot:
        plot_dir = os.path.join(output_dir, 'plots') 
        filedir = os.path.join(plot_dir, '0th_order_y_displacement.png')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir) 
        plt.savefig(filedir, bbox_inches='tight', dpi=300)

    if show_plot:
        plt.show(block=True)

    plt.close() # save memory

    return 


def plot_shifts_cc(exp_time, shifts, ylabel = 'X shift (pixels)',
                   show_plot = False, save_plot = False, 
                   filename = None, output_dir = None):
    """Plots the dispersion or cross-dispersion 
    shifts calculated via cross-correlation

    Args:
        exp_time (np.array): exposure times for each frame.
        x_shifts (np.array): shifts calculated via cross-correlation
        ylabel (str, optional): y-axis label. Defaults to 'X shift (pixels)'.
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
    plt.plot(exp_time, shifts, '-o', color='indianred')
    plt.xlabel('Time Of Exposure (MJD)')
    plt.ylabel(ylabel)
   
    if save_plot:
        plot_dir = os.path.join(output_dir, "plots")
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        plt.savefig(os.path.join(plot_dir, filename),
                    dpi=300,bbox_inches='tight')
        
    if show_plot:
        plt.show(block=True)
    
    plt.close() # save memory

    return 

