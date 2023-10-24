import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.constants import si
from astropy.table import QTable

from specutils import Spectrum1D, SpectralRegion
from specutils.analysis import equivalent_width, snr_derived
from specutils.fitting import fit_generic_continuum, find_lines_derivative, find_lines_threshold
from specutils.manipulation import noise_region_uncertainty


os.path.abspath(os.getcwd())


def read_spectrum(n, filename):
    """
    Function to read a spectrum from a file

    Parameters
    ----------
    n : int
        Number of files to read
    filename : list
        List of files to read

    Returns
    -------
    data : pandas.DataFrame
        Data read from the file
    """
    try:
        data = pd.concat((pd.read_csv('./' + f, sep='\s+', header=None, names=['rx(km/s)', 'ry(Tmb)']) for f in filename), ignore_index=True)
        #for i in range(n):
        #    dat = pd.read_csv(filename[i], sep='\s+', header=None, names=['rx(km/s)', 'ry(Tmb)'])
        #    data = dat if i == 0 else pd.concat([data, dat], ignore_index=True)
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
    

def read_synthetic_spectra(filename):
    """
    Function to read a synthetic spectrum from a file

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
        path = '/Users/teresahuertas/Documents/IAC/Proyectos/Espectros_Sinteticos/GaoLei_Spectra/'
        data = pd.read_csv(path + filename, sep='\s+', header=None, names=['Freq[MHz]', 'Aij'], skiprows=8)
        data.astype({'Freq[MHz]': 'float64', 'Aij': 'float64'}).dtypes
        return data.reset_index(drop=True)
    except FileNotFoundError:
        print(f'File {filename} not found')
        return None
    except Exception as e:
        print(f'Error reading file {filename}: {e}')
        return None


def create_spectrum(data, restfreq, vel, offset=None):
    """
    Function to create a spectrum from a data frame

    Parameters
    ----------
    data : pandas.DataFrame
        Data read from the file
    restfreq : float
        Rest frequency of the spectrum
    vel : float
        Velocity of the source
    offset : float, optional
        Offset of the spectrum
        Default: None
        
    Returns
    -------
    spectrum : specutils.Spectrum1D
        Spectrum created from the data
    """
    # Define equivalence velocity - frequency
    vel_to_freq = [(u.km/u.s, u.MHz, 
                    lambda x: (1-x/si.c.to_value('km/s')) * (restfreq * u.MHz),
                    lambda x: (restfreq * u.MHz -x) / (restfreq * u.MHz) * si.c.to_value('km/s'))]
    
    # Set units
    frequency = (data['rx(km/s)'].values * u.km / u.s).to(u.MHz, equivalencies=vel_to_freq)
    flux = data['ry(Tmb)'].values * u.K

    # Check if restfreq and vel are given with units or not
    if not isinstance(restfreq, u.quantity.Quantity):
        restfreq = restfreq * u.MHz
    if not isinstance(vel, u.quantity.Quantity):
        vel = vel * u.km / u.s
    
    offset = 0.0 * u.MHz if offset is None else offset * u.MHz

    # Create spectrum
    spectrum = Spectrum1D(flux=flux, spectral_axis=frequency+offset, 
                          velocity_convention='radio', 
                          rest_value=restfreq, 
                          radial_velocity=vel)
    
    print('snr: ', snr_derived(spectrum))

    return spectrum


def line_inspector(spectrum, rms, source, line_type=None):
    """
    Function to find lines in a spectrum

    Parameters
    ----------
    spectrum : specutils.Spectrum1D
        Spectrum to find lines
    rms : float
        RMS of the spectrum
    line_type : str, optional
        Type of lines to find (emission or absorption)

    Returns
    -------
    lines : QTable
        Table with the lines found
    """
    # Set noise region
    noise_region = SpectralRegion(min(spectrum.spectral_axis), 
                                  max(spectrum.spectral_axis))
    lines = find_lines_derivative(spectrum, flux_threshold=rms)
    clean_lines = QTable(names=('line_center', 'line_type', 'line_center_index'),
                            dtype=('float64', 'str', 'int64'))

    if line_type == 'emission':
        # Check line width for each emission line
        freq_width = line_freq_width(87317.0)
        for i in range(1, len(lines)):
            if abs(lines[i]['line_center'] - lines[i-1]['line_center']) > freq_width/2:
                print(i, lines[i]['line_center'])
                clean_lines.add_row([lines[i]['line_center'], 
                                     lines[i]['line_type'], 
                                     lines[i]['line_center_index']])
            else:
                print(i, 'No')
            #if equivalent_width(spectrum, regions=SpectralRegion(lines[i]['line_center']-10 * u.MHz, 
            #                                                     lines[i]['line_center']+10 * u.MHz)) < freq_width:
            #    lines.remove_row(i)
        # Save lines to file
        clean_lines.write(f'{source}_emission_lines.txt', format='ascii.ecsv', overwrite=True)
        #print(lines[1]['line_center'])
        #print(lines.info)
        return clean_lines[clean_lines['line_type'] == 'emission'] #lines[lines['line_type'] == 'emission']
    elif line_type == 'absorption':
        # Save lines to file
        lines.write(f'{source}_absorption_lines.txt', format='ascii.ecsv', overwrite=True)
        return lines[lines['line_type'] == 'absorption']
    else:
        # Save lines to file
        lines.write(f'{source}_lines.txt', format='ascii.ecsv', overwrite=True)
        return lines
    

def plot_spectrum(spectrum, lines=None):
    """
    Function to plot a spectrum

    Parameters
    ----------
    spectrum : specutils.Spectrum1D
        Spectrum to plot
    lines : QTable, optional
        Table with the detected lines to plot
    """
    fig, ax = plt.subplots()
    ax.plot(spectrum.spectral_axis, spectrum.flux)
    ax.hlines(0.0, min(spectrum.spectral_axis/u.MHz), 
              max(spectrum.spectral_axis/u.MHz), color='black')
    if lines is not None:
        ax.vlines([lines['line_center']], 0, max(spectrum.flux/u.K), colors='r')
    ax.set_xlabel(spectrum.spectral_axis.unit)
    ax.set_ylabel(spectrum.flux.unit)
    plt.show()

    return None

def plot_synthetic_spectrum(spectrum, lines):
    """
    Function to plot the synthetic data with the observational spectrum

    Parameters
    ----------
    spectrum : specutils.Spectrum1D
        Spectrum to plot
    lines : pandas.DataFrame
        Data with the detected lines to plot

    Returns
    -------
    None
    """
    fig, ax = plt.subplots()
    ax.plot(spectrum.spectral_axis, spectrum.flux)
    ax.hlines(0.0, min(spectrum.spectral_axis/u.MHz), 
              max(spectrum.spectral_axis/u.MHz), color='black')
    ax.vlines(lines['Freq[MHz]'], 0, max(spectrum.flux/u.K), colors='r')
    # Set x-axis limits
    ax.set_xlim(min(spectrum.spectral_axis/u.MHz), 
                max(spectrum.spectral_axis/u.MHz))
    ax.set_xlabel(spectrum.spectral_axis.unit)
    ax.set_ylabel(spectrum.flux.unit)
    plt.show()

    return None


def line_freq_width(restfreq, vel_width=None):
    """
    Function to calculate the line width in units of frequency

    Parameters
    ----------
    restfreq : float
        Rest frequency of the line
    vel_width : float, optional
        Velocity width of the line
        Default: 50 km/s

    Returns
    -------
    freq_width : float
        Line width in units of frequency
    """
    # Check if restfreq is given with units or not
    if not isinstance(restfreq, u.quantity.Quantity):
        restfreq = restfreq * u.MHz

    if vel_width == None:
        vel_width = 50.0 * u.km / u.s
    else:
        # Check if vel_width is given with units or not
        if not isinstance(vel_width, u.quantity.Quantity):
            vel_width = vel_width * u.km / u.s

    freq_width = abs(restfreq * vel_width / (si.c.to_value('km/s') * u.km / u.s))

    return freq_width