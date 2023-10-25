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
        data = pd.read_csv(path + filename, sep='\s+', header=None, names=['Freq[MHz]', 'Aij'], skiprows=8)
        data.astype({'Freq[MHz]': 'float64', 'Aij': 'float64'}).dtypes
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
                ax.text(rrls['Freq[MHz]'][i], -0.02, 
                        rrls['Species'][i], rotation=90, fontsize=12, color='red')
        if ufs is not None:
            for i in range(len(ufs)):
                ax.text(ufs['Freq[MHz]'][i], -0.02, 
                        ufs['Species'][i], rotation=90, fontsize=12, color='green')
        if molecules is not None:
            for i in range(len(molecules)):
                ax.text(molecules['Freq[MHz]'][i], -0.03, 
                        molecules['Species'][i], rotation=90, fontsize=12, color='black')
    else:
        pass

    # Plot horizontal grey line at y=0
    ax.axhline(y=0, color='grey', linestyle='--', linewidth=1.2)

    for key, value, value_index in zip(molec_spectra.keys(), molec_spectra.values(), range(len(molec_spectra))):
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