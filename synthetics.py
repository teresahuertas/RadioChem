import pandas as pd


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