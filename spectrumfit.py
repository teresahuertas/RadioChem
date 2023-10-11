import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.constants import si

from specutils import Spectrum1D, SpectralRegion
from specutils.fitting import fit_generic_continuum, find_lines_derivative, find_lines_threshold
from specutils.manipulation import noise_region_uncertainty


def read_spectrum(filename):
    """
    Function to read a spectrum from a file

    Parameters
    ----------
    filename : str
        Name of the file to read

    Returns
    -------
    data : pandas.DataFrame
        Data read from the file
    """
    try:
        data = pd.read_csv(filename, sep='\s+', header=None, names=['rx(km/s)', 'ry(Tmb)'])
        data.astype({'rx(km/s)': 'float64', 'ry(Tmb)': 'float64'}).dtypes
        #data = data[data['ry(Tmb)'] > -1.0]
        # Remove lines with values out of range [-1, 1]
        data = data[(data['ry(Tmb)'] > -1.0) & (data['ry(Tmb)'] < 1.0)]
        return data.reset_index(drop=True)
    except FileNotFoundError:
        print(f'File {filename} not found')
        return None
    except Exception as e:
        print(f'Error reading file {filename}: {e}')
        return None


def set_units(data, restfreq):
    """
    Function to set the units of the data

    Parameters
    ----------
    data : pandas.DataFrame
        Data read from the file
    restfreq : float
        Rest frequency of the spectrum
        
    Returns
    -------
    data : pandas.DataFrame
        Data with units set
    """
    # Define equivalence velocity - frequency
    vel_to_freq = [(u.km/u.s, u.MHz, 
                    lambda x: (1-x/si.c.to_value('km/s')) * (restfreq * u.MHz),
                    lambda x: (restfreq * u.MHz -x) / (restfreq * u.MHz) * si.c.to_value('km/s'))]
    
    # Set units
    data['rx(km/s)'] = (data['rx(km/s)'].values * u.km / u.s).to(u.MHz, equivalencies=vel_to_freq)
    data['ry(Tmb)'] = data['ry(Tmb)'].values * u.K

    return data