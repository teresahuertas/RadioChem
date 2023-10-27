import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from astropy import constants as const
from astropy import units as u
from specutils import Spectrum1D


def read_synthetic_spectra(path, filename):
    """
    Function to read a synthetic spectrum from a file

    Parameters
    ----------
    path : str
        Path to the file
    filename : str
        Name of the file to read

    Returns
    -------
    data : pandas.DataFrame
        Data read from the file
    """
    try:
        data = pd.read_csv(path + filename, sep='\s+', header=None, names=['Freq[MHz]', 'Aij'], skiprows=8)
        data.astype({'Freq[MHz]': 'float64', 'Aij': 'float64'}).dtypes
        # Remove all rows with data['Aij'] == 0
        data = data[data['Aij'] != 0]
        return data.reset_index(drop=True)
    except FileNotFoundError:
        print(f'File {filename} not found')
        return None
    except Exception as e:
        print(f'Error reading file {filename}: {e}')
        return None
    

def plot_synthetic_spectrum(source, obs_spectra, molec_spectra, rrls=None, ufs=None, molecules=None, names=None):
    """
    Function to plot the observational and synthetic spectra

    Parameters
    ----------
    source : str
        Source name
    obs_spectra : specutils.Spectrum1D
        Observational spectrum
    molec_spectra : dictionary of dataframes
        Dictionary of synthetic spectra, each of them stored in a dataframe
    rrls : dataframe, optional
        Dataframe with the RRLs frequencies
    ufs : dataframe, optional
        Dataframe with the UFs frequencies
    molecules : dataframe, optional
        Dataframe with the molecules frequencies
    names : boolean, optional
        Boolean to indicate whether to plot the names of the species or not

    Returns
    -------
    None
    """
    # Get colormap and create a list of colors
    #cmap = get_cmap('Dark2')
    cmap = get_cmap('tab20')
    #colors = [cmap(i) for i in range(len(molec_spectra))]
    values = [i / len(molec_spectra) for i in range(len(molec_spectra))]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(obs_spectra.spectral_axis, obs_spectra.flux, 
            label=source, color='navy', linewidth=1.1, alpha=0.8)
    
    if rrls is not None:
        ax.vlines(rrls['Freq[MHz]'], 0, max(obs_spectra.flux/obs_spectra.flux.unit),
                    color='red', linewidth=1.2, linestyle='--', label='RRLs')

    if ufs is not None:
        ax.vlines(ufs['Freq[MHz]'], 0, max(obs_spectra.flux/obs_spectra.flux.unit),
                    color='green', linewidth=1.2, linestyle='--', label='UFs')

    if molecules is not None:
        ax.vlines(molecules['Freq[MHz]'], 0, max(obs_spectra.flux/obs_spectra.flux.unit),
                    color='black', linewidth=1.2, linestyle='--', label='Molecules')

    if names is True:
        if rrls is not None:
            for i in range(len(rrls)):
                ax.text(rrls['Freq[MHz]'][i], -0.01, 
                        rrls['Species'][i], rotation=90, fontsize=12, color='red')
        if ufs is not None:
            for i in range(len(ufs)):
                ax.text(ufs['Freq[MHz]'][i], -0.01, 
                        ufs['Species'][i], rotation=90, fontsize=12, color='green')
        if molecules is not None:
            for i in range(len(molecules)):
                ax.text(molecules['Freq[MHz]'][i], -0.01, 
                        molecules['Species'][i], rotation=90, fontsize=12, color='black')
    else:
        pass

    # Plot horizontal grey line at y=0
    ax.axhline(y=0, color='grey', linestyle='--', linewidth=1.2)

    for key, value, value_index in zip(molec_spectra.keys(), 
                                       molec_spectra.values(), 
                                       range(len(molec_spectra))):
        # in molec_spectra.items():
        color = cmap(values[value_index])
        ax.vlines(value['Freq[MHz]'], 0, max(obs_spectra.flux/obs_spectra.flux.unit),
                  color=color, linewidth=1.0, label=key)
    
    # Set axis limits
    ax.set_xlim(min(obs_spectra.spectral_axis/obs_spectra.spectral_axis.unit), 
                max(obs_spectra.spectral_axis/obs_spectra.spectral_axis.unit))
    
    # Set legend outside the plot
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    
    ax.set_xlabel(obs_spectra.spectral_axis.unit)
    ax.set_ylabel(obs_spectra.flux.unit)
    #ax.set_title('Spectrum')
    ax.legend()

    plt.show()


def anttemp(molec_spectra, width, dilution=None, tau=None):
    """
    Function to calculate the antenna temperature of a synthetic spectrum
    
    Parameters
    ----------
    molec_spectra : dictionary of dataframes
        Dictionary of synthetic spectra, each of them stored in a dataframe
    width : float
        Width of the line in velocity units
    dilution : float, optional
        Dilution factor
    tau : float, optional
        Optical depth

    Returns
    -------
    None
    """
    if tau is None:
        corr_tau = 1.0
    else:
        corr_tau = tau / (1 - np.exp(-tau))
    if dilution is None:
        dilution = 1.0
    population = 1.0
    c3 = const.c * const.c * const.c * 1e-3 * 1e-3 * 1e-3

    for key, value in molec_spectra.items():
        # Calculate the antenna temperature
        value['T_A'] = const.h * c3 / (8 * np.pi * const.k_B * value['Freq[MHz]'] * value['Freq[MHz]'] * width) * population * value['Aij'] * dilution / corr_tau
        # Store results in a new column in the dataframe
        molec_spectra[key] = value

    return molec_spectra


def create_synth_spectrum(molec_spectra):
    """
    Function to create a synthetic spectrum from a dictionary of dataframes

    Parameters
    ----------
    molec_spectra : dictionary of dataframes
        Dictionary of synthetic spectra, each of them stored in a dataframe

    Returns
    -------
    synth_spectra : dictionary of specutils.Spectrum1D
        Dictionary of synthetic spectra, each of them stored in a Spectrum1D object
    """
    synth_spectra = {}

    for key, value in molec_spectra.items():
        # Create a Spectrum1D object from the dataframe
        spectrum = Spectrum1D(spectral_axis=value['Freq[MHz]'].values * u.MHz, flux=value['T_A'].values * u.K)
        # Store the Spectrum1D object in the dictionary
        synth_spectra[key] = spectrum

    return synth_spectra