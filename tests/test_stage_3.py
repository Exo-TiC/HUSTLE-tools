import os
import shutil
from urllib import request
import unittest
from astroquery.mast import Observations
from astropy.io import fits
import numpy as np
from scipy.signal import medfilt2d
import xarray as xr

from hustle_tools import stage_0, stage_1, stage_2, stage_3

class TestStage3(unittest.TestCase):
    """ Test hustle_tools stage 3. """

    @classmethod
    def setUpClass(cls):
        # Define local test path, and clear cached data.
        cls.local_data_path = 'test_hustle_tools_data'
        if os.path.exists(cls.local_data_path):
            shutil.rmtree(cls.local_data_path)
        os.mkdir(cls.local_data_path)

        # Create mock 1D spectral dataset.
        np.random.seed(1337) # to keep random draws consistent, apply a seed.
        wavelengths = np.arange(1,11,1) # covers waves 1 to 10
        x = np.array([i for i, z in enumerate(wavelengths)]) # ten columns in this 1D spectrum
        exp_time, orbit_numbers, spec_disp, prof_disp, widths = np.empty((25,)),np.empty((25,)),np.empty((25,)),np.empty((25,10)),np.empty((25,10))
        t0, dt = 0, 5
        ind = 0
        for i in range(5):
            for j in range(5):
                exp_time[ind] = t0
                orbit_numbers[ind] = i+1
                spec_disp[ind] = j
                prof_disp[ind,:] = j
                widths[ind,:] = np.random.normal(1,0.1)
                t0 += dt
                ind += 1
            t0 += 25
        # Fake the y pos.
        y = np.ones((len(exp_time),len(x)))

        # Now build the 1D spectra.
        spec = np.empty((len(exp_time),len(x)))
        spec_err = np.empty_like(spec)
        for i in range(len(x)):
            spec[:,i] = np.random.normal(1,0.05,25)
            spec_err[:,i] = np.array([0.01 for i in spec])
            for j in range(25):
                orbit_n = orbit_numbers[j]
                if orbit_n in (2,3,4):
                    spec[j,i] *= 0.99

        # Compile into an xarray.
        spectra = xr.Dataset(
            data_vars=dict(
                spec = (['exp_time', 'x'], spec),
                spec_err = (['exp_time', 'x'], spec_err),
                trace = (['exp_time', 'x'], y),
                meanstar_disp = (['exp_time', 'x'], y),
                orbit_numbers = (['exp_time',], orbit_numbers),
                spec_disp = (['exp_time'], spec_disp),
                prof_disp = (['exp_time', 'x'], prof_disp),
                widths = (['exp_time', 'x'], widths),
                ),
            coords=dict(
                wave=(['x'], wavelengths),
                trace_x = x,
                exp_time = exp_time,
                ),
            )

        # Download lightweight test data. First two exposures of first two
        # orbits of W31, proposal_id=17183, visit=16.
        cls.visit_number = "16"
        cls.mast_data_files = [
            "iexr16ljq_flt.fits",  # Direct images.
            "iexr16ljq_jit.fits",
            "iexr16ljq_spt.fits",
            "iexr16lkq_flt.fits",  # Spec images.
            "iexr16lkq_jit.fits",
            "iexr16lkq_spt.fits",
            "iexr16llq_flt.fits",
            "iexr16llq_jit.fits",
            "iexr16llq_spt.fits",
            "iexr16luq_flt.fits",
            "iexr16luq_jit.fits",
            "iexr16luq_spt.fits",
            "iexr16lvq_flt.fits",
            "iexr16lvq_jit.fits",
            "iexr16lvq_spt.fits"]
        for mdf in cls.mast_data_files:
            Observations.download_file(
                "mast:HST/product/{}".format(mdf),
                local_path=os.path.join(cls.local_data_path, mdf))
        stage_0.collect_and_move_files(
            cls.visit_number, cls.local_data_path, cls.local_data_path)
        # Save spectra xarray out.
        spectra.to_netcdf(os.path.join(cls.local_data_path, 'specs_+1.nc'))

        cls.spectra = spectra

        print("Test data is ready to go!")

    @classmethod
    def tearDownClass(cls):
        # Tidy up.
        if os.path.exists(cls.local_data_path):
            shutil.rmtree(cls.local_data_path)
    
    def test_a_load_data(self):
        """ Load .nc files as xarray."""
        spectra = stage_3.load_data_S3(self.local_data_path, order="+1")
        self.assertEqual(spectra.spec.shape[0], 25)

    def test_b_bin_light_curves(self):
        """ Bin 1D spectra into light curves. """
        # First, test binning by wavelengths. We will get five light curves from this.
        light_curves = stage_3.bin_light_curves(self.spectra, 
                                                "+1", 
                                                bin_method = "wavelengths",
                                                bins = np.array([0,2,4,6,8,10]),
                                                ncol = 2,
                                                normalize = True,
                                                norm_lim = len([i for i in self.spectra.orbit_numbers.data if i == 1]), 
                                                rem_exp = None)
        
        # Check that the binned results are what we expect.
        self.assertEqual(light_curves.spec_lc.shape[0],5,
                         msg='Wavelength binning failed!')
        self.assertEqual(light_curves.spec_lc.shape[1],25)

        # Now bin by columns. We again expect to get five light curves here.
        light_curves = stage_3.bin_light_curves(self.spectra, 
                                                "+1", 
                                                bin_method = "columns",
                                                bins = np.array([0,2,4,6,8,10]),
                                                ncol = 2,
                                                normalize = True,
                                                norm_lim = len([i for i in self.spectra.orbit_numbers.data if i == 1]), 
                                                rem_exp = None)
        
        # Check that the binned results are what we expect.
        self.assertEqual(light_curves.spec_lc.shape[0],5,
                         msg='Column binning failed!')
        self.assertEqual(light_curves.spec_lc.shape[1],25)

    def test_c_read_vectors(self):
        """ Read *jit.fits files for jitter information. """
        # WIP!
        state_vectors = stage_3.get_state_vectors(self.spectra,
                                                  data_dir = self.local_data_path,
                                                  method = 'jit_dec',
                                                  include_jitter = False,
                                                  plot = False)

if __name__ == '__main__':
    unittest.main()