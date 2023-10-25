import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
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
        #path = '/Users/teresahuertas/Documents/IAC/Proyectos/Espectros_Sinteticos/GaoLei_Spectra/'
        data = pd.read_csv(path + filename, sep='\s+', header=None, names=['Freq[MHz]', 'Aij'], skiprows=8)
        data.astype({'Freq[MHz]': 'float64', 'Aij': 'float64'}).dtypes
        return data.reset_index(drop=True)
    except FileNotFoundError:
        print(f'File {filename} not found')
        return None
    except Exception as e:
        print(f'Error reading file {filename}: {e}')
        return None
    

def plot_synthetic_spectrum(source, obs_spectra, molec_spectra, rrls=None, ufs=None):
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
        ax.vlines(rrls['Freq[MHz]'], 0, max(obs_spectra.flux/u.K),
                    color='red', linewidth=1.0, label='RRLs')
    if ufs is not None:
        ax.vlines(ufs['Freq[MHz]'], 0, max(obs_spectra.flux/u.K),
                    color='green', linewidth=1.0, label='UFs')
    # Plot horizontal grey line at y=0
    ax.axhline(y=0, color='grey', linestyle='--', linewidth=1.2)
    for key, value, value_index in zip(molec_spectra.keys(), molec_spectra.values(), range(len(molec_spectra))):
        # in molec_spectra.items():
        color = cmap(values[value_index])
        ax.vlines(value['Freq[MHz]'], 0, max(obs_spectra.flux/u.K),
                  color=color, linewidth=1.0, label=key)
    
    # Set axis limits
    ax.set_xlim(min(obs_spectra.spectral_axis/u.MHz), 
                max(obs_spectra.spectral_axis/u.MHz))
    
    # Set legend outside the plot
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    
    ax.set_xlabel(obs_spectra.spectral_axis.unit)
    ax.set_ylabel(obs_spectra.flux.unit)
    #ax.set_title('Spectrum')
    ax.legend()

    plt.show()