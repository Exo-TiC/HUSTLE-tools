import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


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
    
    raw_wlc = np.sum(spec, axis=1)

    plt.figure(figsize = (10, 7))
    plt.plot(times, raw_wlc, 'o', color='indianred', markeredgecolor='black')
    plt.xlabel('Time of exposure')
    plt.ylabel('Counts')

      
    if save_plot:
        plot_dir = os.path.join(output_dir, 'plots') 
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir) 
        filedir = os.path.join(plot_dir, f'{filename}.png')
        plt.savefig(filedir,dpi=300,bbox_inches='tight')

    if show_plot:
        plt.show(block=True)

    plt.close() # save memory


    return 0


def plot_raw_spectrallightcurves(times, spec, order="+1",
                    show_plot = False, save_plot = False,
                    filename = None, output_dir = None):
    

    for i, lc in enumerate(np.transpose(spec)):
        raw_lc = lc/np.median(lc[:30])

        plt.figure(figsize = (10, 7))
        plt.plot(times, raw_lc, 'o', color='indianred', markeredgecolor='black')
        plt.xlabel('Time of exposure')
        plt.ylabel('Counts')

      
        if save_plot:
            plot_dir = os.path.join(output_dir, 'plots') 
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir) 
            filedir = os.path.join(plot_dir, f'{filename}.png')
            plt.savefig(filedir,dpi=300,bbox_inches='tight')

        if show_plot:
            plt.show(block=True)

        plt.close() # save memory


    return 0



def plot_aperture_lightcurves(obs, tested_hws, wlcs,  
                               show_plot = False, save_plot = False,
                                filename = None, output_dir = None):

  
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
    
    if save_plot:
        plot_dir = os.path.join(output_dir,'plots')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        plt.savefig(os.path.join(plot_dir, filename),
                    dpi=300,bbox_inches='tight')
    if show_plot:
        plt.show(block=True)
    plt.close()


    return 0