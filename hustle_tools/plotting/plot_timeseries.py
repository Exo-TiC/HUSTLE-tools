import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.cm as cm
import matplotlib.pylab as pl


#define plotting parameters
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('axes', labelsize=14)
plt.rc('legend',**{'fontsize':11})


def plot_flags_per_time(series_x, series_y, style='line',
                        line_data = None, scatter_data = None, 
                        title = None, xlabel=None, ylabel=None,
                        xmin = 0, xmax = 1, ymin = 0, ymax = 1e4,
                        mark_size = 30, line_style='-',
                        show_plot = False, save_plot = False, 
                        filename = None, output_dir = None):
    """Function to plot number of flagged pixels vs time.

    Args:
        series_x (array-like): series of x coordinates for plotting.
        series_y (array-like): series of y coordinates for plotting.
        style (str): options are 'line' or 'scatter'. Defaults to 'line'.
        line_data (array-like, optional): x, y values defining lines
        to overplot on top of series_x, series_y. Defaults to None.
        scatter_data (array-like, optional): x, y values defining scatter
        points to overplot on top of series_x, series_y. Defaults to None.
        title (str, optional): title for the plot. Defaults to None.
        xlabel (str, optional): x axis label. Defaults to None.
        ylabel (str, optional): y axis label. Defaults to None.
        xmin (float, optional): x axis lower limit. Defaults to 0.
        xmax (float, optional): x axis upper limit. Defaults to 1.
        ymin (float, optional): y axis lower limit. Defaults to 0.
        ymax (float, optional): y axis upper limit. Defaults to 1e4.
        mark_size (float, optional): size of scatter points. Defaults to 30.
        line_style (str, optional): mpl style of line. Defaults to '-'.
        show_plot (bool, optional): whether to interrupt execution to show the
        user the plot. Defaults to False.
        save_plot (bool, optional): whether to save this plot to a file.
        Defaults to False.
        filename (list of str, optional): names to give each output file,
        if saving. Defaults to None.
        output_dir (str, optional): where to save the files to. Defaults to None.
    """
    
    for i, (x, y) in enumerate(zip(series_x, series_y)):

        plt.figure(figsize = (20, 4))
        if style == 'scatter':
            plt.scatter(x, y, color='k', s=mark_size)
        elif style == 'line':
            plt.plot(x, y, color='k', ls=line_style)
        if xlabel:
            plt.xlabel(xlabel[i])
        if ylabel:
            plt.ylabel(ylabel[i])
        plt.colorbar()

        if xmin or xmax:
            plt.xlim(xmin, xmax)

        if ymin or ymax:
            plt.ylim(ymin, ymax)

        if line_data:
            for j, line in enumerate(line_data):
                plt.plot(line[0], line[1])

        if scatter_data: 
            plt.scatter(scatter_data[0], scatter_data[1], s = mark_size, color = 'r', marker = '+')

        if title:
            plt.title(title)
        
        if save_plot:
            plot_dir = os.path.join(output_dir, 'plots') 
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir) 
            filedir = os.path.join(plot_dir, f'{filename[i]}.png')
            plt.savefig(filedir, bbox_inches = 'tight', dpi = 300)
        
        if show_plot:
            plt.show(block=True)

        plt.close() # save memory
    
    return


def plot_raw_whitelightcurve(times, spec, order="+1",
                             show_plot = False, save_plot = False,
                             filename = None, output_dir = None):
    """Plots the uncorrected broad-band light curve for this order, as a
    diagnostic of your cleaning process.

    Args:
        times (np.array): mid-exposure time of each frame.
        spec (np.array): 1D extracted spectra.
        order (str, optional): which order we are plotting, for plot title.
        Defaults to "+1".
        show_plot (bool, optional): whether to interrupt execution to show the
        user the plot. Defaults to False.
        save_plot (bool, optional): whether to save this plot to a file.
        Defaults to False.
        filename (str, optional): name to give this file, if saving.
        Defaults to None.
        output_dir (str, optional): where to save the file, if saving.
        Defaults to None.
    """

    # define order colors
    colors = {"+1":'indianred',"-1":'dodgerblue',
              "+2":'orangered',"-2":'royalblue',
              "+3":'darkorange',"-3":'blue',
              "+4":'orange',"-4":'deepskyblue'}

    raw_wlc = np.sum(spec, axis=1)

    plt.figure(figsize = (10, 7))
    plt.plot(times, raw_wlc, 'o', color=colors[order], markeredgecolor='black')
    plt.xlabel('Time of exposure')
    plt.ylabel('Counts')
    plt.title("Raw broad-band light curve, order {}".format(order))

      
    if save_plot:
        plot_dir = os.path.join(output_dir, 'plots') 
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir) 
        filedir = os.path.join(plot_dir, f'{filename}.png')
        plt.savefig(filedir,dpi=300,bbox_inches='tight')

    if show_plot:
        plt.show(block=True)

    plt.close() # save memory

    return 


def plot_raw_spectrallightcurves(times, spec, order="+1",
                                 show_plot = False, save_plot = False,
                                 filename = None, output_dir = None):
    """Plots the uncorrected spectrally-binned light curves for this order, as
    diagnostics of your cleaning process.

    Args:
        times (np.array): mid-exposure time of each frame.
        spec (np.array): 1D extracted spectra.
        order (str, optional): which order we are plotting, for plot title.
        Defaults to "+1".
        show_plot (bool, optional): whether to interrupt execution to show the
        user the plot. Defaults to False.
        save_plot (bool, optional): whether to save this plot to a file.
        Defaults to False.
        filename (str, optional): name to give this file, if saving.
        Defaults to None.
        output_dir (str, optional): where to save the file, if saving.
        Defaults to None.
    """

    # define order colors
    colors = {"+1":'indianred',"-1":'dodgerblue',
              "+2":'orangered',"-2":'royalblue',
              "+3":'darkorange',"-3":'blue',
              "+4":'orange',"-4":'deepskyblue'}    

    for i, lc in enumerate(np.transpose(spec)):
        n_oot = int(0.20*lc.shape[0]) # typically, first 20% of data is the first orbit, which is oot/ooe
        raw_lc = lc/np.median(lc[:n_oot])

        plt.figure(figsize = (10, 7))
        plt.plot(times, raw_lc, 'o', color=colors[order], markeredgecolor='black')
        plt.xlabel('Time of exposure')
        plt.ylabel('Counts')
        plt.title("{}th column's spectral light curve, order {}".format(i,order))

      
        if save_plot:
            plot_dir = os.path.join(output_dir, 'plots') 
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir) 
            filedir = os.path.join(plot_dir, f'{filename}_lc{i}.png')
            plt.savefig(filedir,dpi=300,bbox_inches='tight')

        if show_plot:
            plt.show(block=True)

        plt.close() # save memory

    return 


def plot_aperture_lightcurves(obs, tested_hws, wlcs,  
                              show_plot = False, save_plot = False,
                              filename = None, output_dir = None):
    """Plot each extracted broad-band light curve per halfwidth,
    to show which halfwidth produced the nicest light curve.

    Args:
        obs (xarray): just need the .exp_time from this.
        tested_hws (array-like): int, the halfwidths of extraction
        that we tested.
        wlcs (array-like): each light curve extracted.
        show_plot (bool, optional): whether to interrupt execution to show the
        user the plot. Defaults to False.
        save_plot (bool, optional): whether to save this plot to a file.
        Defaults to False.
        filename (str, optional): name to give this file, if saving.
        Defaults to None.
        output_dir (str, optional): where to save the file, if saving.
        Defaults to None.
    """

    # colormap
    cmap = cm.get_cmap('viridis')
    cs = cmap(np.linspace(0,1,len(tested_hws)))

    # offsets
    #offsets = np.arange(0, len(tested_hws))*0.001

    plt.figure(figsize=(10, 7))
    for wlc, hw, c in zip(wlcs, tested_hws, cs):
        if (hw == tested_hws[0] or hw == tested_hws[-1]):
            plt.scatter(obs.exp_time, wlc, color=c, label=hw, alpha=0.75)
        else:
            plt.scatter(obs.exp_time, wlc, color=c, alpha=0.75)
    plt.legend(loc='upper left', ncols=2)
    plt.xlabel('Time of exposure')
    plt.ylabel('Counts')
    plt.title("Light curve for each tested halfwidth")
    
    if save_plot:
        plot_dir = os.path.join(output_dir,'plots')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        plt.savefig(os.path.join(plot_dir, f'{filename}.png'),
                    dpi=300,bbox_inches='tight')
        
    if show_plot:
        plt.show(block=True)

    plt.close() # save memory

    return


def plot_raw_binned_spectrallightcurves(light_curves, order, show_plot = False, save_plot = False,
                                        filename = None, output_dir = None):
    """Plots a gif of the raw binned spectral light curves.

    Args:
        light_curves (_type_): _description_
        order (_type_): _description_
        show_plot (bool, optional): _description_. Defaults to False.
        save_plot (bool, optional): _description_. Defaults to False.
        filename (_type_, optional): _description_. Defaults to None.
        output_dir (_type_, optional): _description_. Defaults to None.
    """
    
    # get needed data
    times = light_curves.exp_time.data
    spec_lcs = light_curves.spec_lc.data
    wave_cents = light_curves.wave_cents.data

    # define colors
    colors = pl.cm.jet(np.linspace(0, 1, len(wave_cents)))

    # create animation of light curves
    fig,ax = plt.subplots(figsize=(10, 7))
    
    # plot first light curve to get things started
    spec_line, = ax.plot(times,spec_lcs[0],'o',ls='--',color = colors[0], markeredgecolor='black',
                        label=r"{} order, {:.0F} $\AA$".format(order,wave_cents[0]))
    spec_line.set_color(colors[0])
    ax.axhline(y=1,color='k',ls=':')
    leg = ax.legend(loc='upper right')
    #ax.set_xlim(2000,8000)
    #ax.set_ylim(0.99, 1.01)
    ax.set_xlabel('Time of exposure')
    ax.set_ylabel('Relative flux')
    ax.set_title("{} order spectral light curve []".format(order))

    # initialize 
    def init():
        spec_line.set_data([times,spec_lcs[0]])
        spec_line.set_color(colors[0])
        leg.get_texts()[0].set_text(r"{} order, {:.0F} $\AA$".format(order,wave_cents[0]))

        return spec_line,

    # define animation function
    def animation_func(i):
        # update line data
        spec_line.set_data([times,spec_lcs[i]])
        spec_line.set_color(colors[i])
        leg.get_texts()[0].set_text(r"{} order, {:.0F} $\AA$".format(order,wave_cents[i]))

        return spec_line,
        
    # create and plot animation
    animation = FuncAnimation(fig, animation_func, init_func = init,
                              frames = len(spec_lcs), interval = 20)
    plt.tight_layout()

    # save animation
    if save_plot:
        plot_dir = os.path.join(output_dir, 'plots') 

        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

        animation.save(os.path.join(plot_dir, f'{filename}.gif'), writer = 'ffmpeg', fps = 5)

    if show_plot:
        plt.show(block = True)

    plt.close() # save memory

    return 


def plot_waterfall(light_curves, order, show_plot=False, save_plot=False,
                   filename=None, output_dir = None):
    """_summary_

    Args:
        light_curves (_type_): _description_
        order (_type_): _description_
        show_plot (bool, optional): _description_. Defaults to False.
        save_plot (bool, optional): _description_. Defaults to False.
        filename (_type_, optional): _description_. Defaults to None.
        output_dir (_type_, optional): _description_. Defaults to None.
    """

    # import data
    exp_time = light_curves.exp_time.data
    wave_cents = light_curves.wave_cents.data
    wave_edges = light_curves.wave_edges.data
    lcs_raw = light_curves.spec_lc.data
    
    # define colors and offsets
    colors = pl.cm.jet(np.linspace(0, 1, len(wave_cents)))
    offset = -0.02 # find a way to calculate this automatically

    # get index cut
    half_ind = int(np.shape(lcs_raw)[0]/2)
  
    # create figure
    plt.figure(figsize = (8, 35))
    plt.title(f'Raw binned light curves order {order}')
    ax1 = plt.subplot2grid((1, 2), (0, 0))

    for i, lc in enumerate(lcs_raw[:half_ind]):

        ax1.plot(exp_time, lc/np.mean(lc) + i*offset, 
                 '--', color = colors[i], markersize=3, linewidth=1)
        
        ax1.plot(exp_time, lc/np.mean(lc) + i*offset, 
                 'o', color = colors[i], markersize=3, markeredgecolor = 'black', markeredgewidth = 0.5)
        
        ax1.text(exp_time[0], 1 + i*offset - offset/2, 
                 '[{:.0f}, {:.0f}] $\AA$'.format(wave_edges[i], wave_edges[i+1]), 
                 color = colors[i], 
                 fontsize = 8, 
                 fontweight='bold')

         
    ax1.set_ylabel('Normalized Flux')
    ax1.set_xlabel('Time from Mid-transit (days)')
    #ax1.set_ylim(ylims)

    ax2 = plt.subplot2grid((1, 2), (0, 1), sharey = None, sharex = ax1)

    for i, lc in enumerate(lcs_raw[half_ind:]):

        ax2.plot(exp_time, lc/np.mean(lc) + i*offset, 
                 '--', color = colors[half_ind + i], markersize=3, linewidth=1)
        
        ax2.plot(exp_time, lc/np.mean(lc) + i*offset, 
                 'o', color = colors[half_ind + i], markersize=3, markeredgecolor = 'black', markeredgewidth = 0.5)
        
        ax2.text(exp_time[0], 1 + i*offset - offset/2, 
                 '[{:.0f}, {:.0f}] $\AA$'.format(wave_edges[half_ind + i], wave_edges[half_ind + i + 1]), 
                 color = colors[half_ind + i], 
                 fontsize = 8, 
                 fontweight='bold')
     
    plt.subplots_adjust(wspace=0.05, hspace=0.1)

    ax2.set_xlabel('Time from Mid-transit (days)')
    ax2.set_yticks([])


    if save_plot:
        plot_dir = os.path.join(output_dir, 'plots') 
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir) 
        filedir = os.path.join(plot_dir, f'{filename}_order{order}.png')
        plt.savefig(filedir, dpi=300, bbox_inches='tight')

    if show_plot:
        plt.show(block=True)

    plt.close() # save memory


    return 

