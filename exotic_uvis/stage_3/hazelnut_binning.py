import numpy as np
import matplotlib.pyplot as plt
#from scipy.stats import linregress
from scipy.optimize import least_squares

def do(electrons, wavelength_solutions, wavbins):
    '''
    Converts 1D spectra into binned spectroscopic curves.
    
    :param electrons: 3D array. The electrons(#Order,t,lambda) array.
    :param wavelength_solutions: lst of lsts of float. The wavs(wavelength) array in Angstroms for each order.
    :param wavbins: list of floats. The edges defining each spectroscopic light curve. The ith bin will count pixels that have wavelength solution wavbins[i] <= wav < wavbins[i+1].
    :return: curves(Order#, t, lambda) array, central_lam list. The spectroscopic light curves object and corresponding central wavelengths in Angstroms.
    '''

    # Initialize variables.
    curves = np.empty((len(wavelength_solutions),np.shape(electrons)[1],len(wavbins)-1))
    WLC = np.empty((len(wavelength_solutions),np.shape(electrons)[1]))
    #central_lams = []
    central_lams = np.empty((len(wavelength_solutions),len(wavbins)-1))
    
    # Iterate over orders.
    for j, wavs in enumerate(wavelength_solutions):
        # Iterate over time in that order.
        for i in range(np.shape(electrons)[1]):
            # Grab jth order spectrum at time i.
            flx = electrons[j,i,:]
            if len(flx) != len(wavs):
                # If the flx object is too long, truncate it to have the right length as the wavs object.
                # The right edge gets cut off, and the right edge is 0'd out if wavs was not long enough so it's fine.
                flx = flx[:len(wavs)]
            for k in range(len(wavbins)-1):
                # Iterate over central wavelength and get flux element for discrete wavelength bin.
                #central_lams.append((wavbins[k]+wavbins[k+1])/2)
                ok = (wavs>wavbins[k]) & (wavs<wavbins[k+1])
                curves[j,i,k] = np.sum(flx[ok])
            # The WLC element is just the sum of the entire 1D spectrum at that time in that order.
            WLC[j,i] = np.sum(curves[j,i,:])
        for k in range(len(wavbins)-1):
            # Build this order's central wavelengths object.
            central_lams[j,k] = (wavbins[k]+wavbins[k+1])/2
    
    # While the WLC is not median-normalized, this is our one shot to get at the photon noise limit.
    for j, wavs in enumerate(wavelength_solutions):
        med = np.median(WLC[j,:])
        print("Median total counts in the {}th WLC: {}".format(j, med))
        print("Photon noise limit assuming you have removed all non-signal sources: {} ppm".format(10**6 / (med**0.5)))
    print("Extracted a white light curve and {} spectroscopic light curves.".format(len(wavbins)-1))
    return WLC, curves, central_lams

def est_errs(time, flx, kick_outliers=True):
    #result = linregress(x,flx)
    #fit = result.slope*x + result.intercept
    try:
        result = least_squares(residuals_,
                            np.array([1,1,0,1]),
                            args=(time,flx))
    except ValueError:
        plt.scatter(time,flx)
        plt.title("somethin wrong with you")
        plt.show()
        plt.close()
        print(1/0)
    '''
    a,b,c,d = result.x
    plt.scatter(x,flx)
    plt.plot(x,rampslope(x,a,b,c,d))
    plt.show()
    plt.close()
    '''
    res = residuals_(result.x,time,flx,kick_outliers)
    return res

def residuals_(fit,x,flx,kick_outliers=False):
    rs = rampslope(x,fit[0],fit[1],fit[2],fit[3])
    residuals = flx-rs
    if kick_outliers:
        # Remove outlier points.
        res_mean = np.mean(residuals)
        res_sig = np.std(residuals)
        outliers = np.where(np.abs(residuals-res_mean) > 3*res_sig)[0]
        residuals = np.delete(residuals, outliers)
    return residuals

def rampslope(x,a,b,c,d):
    return a*np.exp(b*x) + c*x + d

import numpy as np

def mednorm(curves, indices=None):
    '''
    Median-normalizes the given light curves.
    
    :param curve: 2D or 3D array. If 2D, is the WLC. If 3D, is the spec curves object.
    :param indices: lst of int. Indices of out-of-transit or in-eclipse flux points. Used to normalize, if not None.
    :return: curve 2D/3D array that has been median normalized.
    '''
    # Get dims to iterate over.
    orders  = range(np.shape(curves)[0])
    is_WLC = False
    try:
        lambdas = range(np.shape(curves)[2])
    except:
        # It's the WLC.
        is_WLC = True
    
    for o in orders:
        if is_WLC:
            a = curves[o,:]
            a = a[a != 0]
            if indices:
                # Use median of specified indices. Will need to rebuild a real quick.
                a_t = []
                for ind in indices:
                    for i in a[ind[0]:ind[1]]:
                        a_t.append(i)
                med = np.median(np.array(a_t))
            else:
                # Use median of entire event.
                med = np.median(a)
            curves[o,:] = curves[o,:]/med
        else:
            for l in lambdas:
                a = curves[o,:,l]
                a = a[a != 0]
                med = np.median(a)
                curves[o,:,l] = curves[o,:,l]/med
    return curves

import os

def write_spec_curve(curves, timestamps, errors, orders, wavbins, wavelength_solutions, outdir):
    '''
    Write every spec curve in curves object to the outdir.

    :param curves: 3D array. Array of [o,t,l] fluxes to write out.
    :param timestamps: 1D array. Timestamps in MJD for flux points.
    :param errors: 3D array. Array of [o,t,l] errors to write out.
    :param orders: lst of str. Orders that are in the array (e.g. "p1", "m1", "p2", "m2").
    :param wavbins: lst of float. Edges of bins used to build light curves.
    :param wavelength_solutions: lst of lsts. Wavelength solutions for each order.
    :param outdir: str. Where to write .txt files to.
    :return: outfiles written to outdir, with subdirectories for each order.
    '''
    for i, order in enumerate(orders):
        # Iterate over order.
        orderdir = outdir + "_{}".format(order)
        for j in range(curves.shape[2]):
            # Iterate over wavelength in that order.
            outfile = "slc_%.0f_%.0f" % (wavbins[j],wavbins[j+1])
            write_one_curve(timestamps,curves[i,:,j],errors[i,:,j],outfile,orderdir)
            outfile = "wvs_%.0f_%.0f" % (wavbins[j],wavbins[j+1])
            write_one_waves([wavbins[j],wavbins[j+1]], wavelength_solutions[i], outfile, orderdir)

def write_wl_curve(curves, timestamps, errors, orders, wavbins, wavelength_solutions, outdir):
    '''
    Write every white light curve in curves object to the outdir.

    :param curves: 2D array. Array of [o,t] fluxes to write out.
    :param timestamps: 1D array. Timestamps in MJD for flux points.
    :param errors: 2D array. Array of [o,t] uncertainties to write out.
    :param orders: lst of str. Orders that are in the array (e.g. "p1", "m1", "p2", "m2").
    :param wavbins: lst of float. Edges of bins used to build white light curve.
    :param wavelength_solutions: lst of lsts. Wavelength solutions for each order.
    :param outdir: str. Where to write .txt files to.
    :return: outfiles written to outdir, with subdirectories for each order.
    '''
    for i, order in enumerate(orders):
        # Iterate over order.
        orderdir = outdir + "_{}".format(order)
        outfile = "wlc"
        write_one_curve(timestamps,curves[i,:],errors[i,:],outfile,orderdir)
        outfile = "wvs_wlc"
        write_one_waves([min(wavbins),max(wavbins)], wavelength_solutions[i], outfile, orderdir)

def write_one_curve(timestamps, curve, errors, outfile, outdir, failed=False):
    '''
    Writes one light curve to .txt file for future reading.

    :param timestamps: 1D array. Timestamps in MJD for flux points.
    :param curve: 1D array. Normalized flux values.
    :param errors: 1D array. Normalized flux errors.
    :param outfile: str. Name of the file to save to.
    :param outdir: str. Where to write the outfile to.
    :return: outfile of timestamps and curve written to outdir/outfile.
    '''
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    with open('{}/{}.txt'.format(outdir,outfile), mode='w') as file:
        if failed:
            file.write("#failed!\n")
        file.write('#time[MJD]     flux[normalized]     uncertainty[normalized]     shift[pixels]\n')
        for t, flx, err in zip(timestamps, curve, errors):
            if failed:
                flx = 1
            file.write('{}     {}     {}     {}\n'.format(t,flx,err,0))

def write_one_waves(wavbins, wavelength_solution, outfile, outdir):
    '''
    Writes one set of wavelengths used to build one light curve to .txt file for future reading.

    :param wavbins: 1D array. The min and max wavelength used to build a light curve.
    :param wavelength_solution: 1D array. The wavelength solution for the x axis.
    :param outfile: str. Name of the file to save to.
    :param outdir: str. Where to write the outfile to.
    :return: outfile of wavelengths written to outdir/outfile.
    '''
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    with open('{}/{}.txt'.format(outdir,outfile), mode='w') as file:
        file.write("#wavelengths[AA]\n")
        for wav in wavelength_solution:
            if (wav < wavbins[1] and wav > wavbins[0]):
                file.write('{}\n'.format(wav))