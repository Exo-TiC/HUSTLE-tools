import numpy as np
from astropy.io import fits
from wfc3tools import sub2full
import matplotlib.pyplot as plt
import xarray as xr
from tqdm import tqdm
import os


def load_data_S3(data_dir, order="+1", verbose = 0):
     
    """

    Function to read in the outputs from stage 2

    """

    specs = xr.open_dataset(os.path.join(data_dir, f'specs_{order}.nc')) 

    return specs

def save_data_S3(light_curves,
                output_dir = None, order='+1',
                filename = 'light_curves'):
    """Function to create and save xarray containing the information extracted
    from stage 3.

    Args:
        obs (_type_): _description_
        specs (_type_): _description_
        specs_err (_type_): _description_
        trace_x (_type_): _description_
        trace_y (_type_): _description_
        wavelengths (_type_): _description_
        orders (tuple, optional): _description_. Defaults to ("+1", "-1").
        output_dir (_type_, optional): _description_. Defaults to None.
        filename (str, optional): _description_. Defaults to 'specs'.

    Returns:
        _type_: _description_
    """

    
    # Save results in Stage 2 folder 
    light_curves.to_netcdf(os.path.join(output_dir, f'{filename}_{order}.nc'))

    return 0
