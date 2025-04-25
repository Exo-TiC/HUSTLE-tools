import os

import xarray as xr


def load_data_S3(data_dir, order="+1"):
    """Reads in the Stage 2 outputs for the requested order.

    Args:
        data_dir (str): filepath pointing to a Stage 2 specs_{order}.nc file.
        order (str, optional): the order you want to process. Defaults to "+1".

    Returns:
        xarray: the loaded Stage 2 1D spectral file for the given order.
    """
    specs = xr.open_dataset(os.path.join(data_dir, f'specs_{order}.nc')) 

    return specs


def save_data_S3(light_curves,
                 output_dir = None, order='+1',
                 filename = 'light_curves'):
    """Function to create and save xarray containing the information extracted
    from stage 3.

    Args:
        light_curves (xarray): xarray containing all binned white light and
        spectroscopic light curves.
        output_dir (str, optional): path to folder where the light curves
        will be saved. Defaults to None.
        order (str, optional): which order this is, for file naming.
        Defaults for '+1'.
        filename (str, optional): name the file will receive, sans order tag.
        Defaults to 'light_curves'.
    """

    # Save results in Stage 3 folder 
    light_curves.to_netcdf(os.path.join(output_dir, f'{filename}_{order}.nc'))

    return 
