import pandas as pd
import os


os.path.abspath(os.getcwd())


class Catalogue:
    """
    Class to read, classify and analyse the obseved lines in PNe
    """

    # Paths to catalogues defined from the current directory
    cat_path = 'RadioChem/Source_Catalogues/'
    results_path = 'RadioChem/Results/'

    def __init__(self, name=None):
        if not name:
            print('Line catalogue created with any source specified')
        elif name == None:
            self.catalogue = None
        else:
            self.catalogue = self.read_catalogue_file(name)
            '''if name == 'rrls':
                #self.catalogue = self.catalogue
            else:
                #self.rrls = 
                #self.molecules =
                #self.uf =
'''
        return
    

    def read_catalogue_file(self, name):
        """
        Function to read the catalogue file and return a pandas dataframe
        
        Parameters
        ----------
        name : str
            Name of the catalogue file to read

        Returns
        -------
        catalogue : pandas dataframe
            Dataframe with the catalogue information
        """
        try:
            file_format = [(0, 1), (1, 4), (4, 19), (19, 29), 
                           (29, 37), (37, 45), (45, 50), (50, 60), 
                           (60, 75), (75, 79), (79, 94), (94, 100)]
            data = pd.read_fwf(name, colspecs=file_format, header=None)
            data = data[5:]
            data.drop([1, 4, 5, 6, 7, 9], axis=1, inplace=True)
            data.rename(columns={0: 'Status',
                                 2: 'Species',
                                 3: 'Freq[MHz]',
                                 8: 'Upper',
                                 10: 'Lower',
                                 11: 'Origin'}, inplace=True)
            print(name + ' data read successfully')
            return data
        except FileNotFoundError:
            print('File not found')
            return None
        except Exception as e:
            print('Error reading {e}')
            return None