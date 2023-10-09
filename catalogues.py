import pandas as pd
import os


os.path.abspath(os.getcwd())


class Catalogue:
    """
    Class to read, classify and analyse the obseved lines in PNe
    """

    # Paths to catalogues defined from the current directory
    cat_path = './Source_Catalogues/'
    results_path = 'RadioChem/Results/'

    def __init__(self, name=None):
        if not name:
            print('Line catalogue created with any source specified')
        elif name == None:
            self.catalogue = None
        else:
            self.catalogue = self.read_catalogue_file(name)
            if name == 'rrls':
                #self.catalogue = self.catalogue
            '''else:
                #self.rrls = 
                #self.molecules =
                #self.uf =
'''
        return
    

    def read_catalogue_file(self, name, cat_path=cat_path):
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
            data = pd.read_fwf(cat_path + name + '.my-lines.list', colspecs=file_format, header=None)
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
    

    def sort_greek(self, data):
        """
        Function to sort lines by greek letter

        Parameters
        ----------
        data : pandas dataframe
            Dataframe with the catalogue information

        Returns
        -------
        final_data : pandas dataframe
            Dataframe with the catalogue information sorted by greek letter
        """
        ga = data[data['Species'].str.endswith('a')]
        gb = data[data['Species'].str.endswith('b')]
        gg = data[data['Species'].str.endswith('g')]
        gd = data[data['Species'].str.endswith('d')]
        ge = data[data['Species'].str.endswith('e')]
        gz = data[data['Species'].str.endswith('z')]
        gh = data[data['Species'].str.endswith('h')]
        gq = data[data['Species'].str.endswith('q')]
        gi = data[data['Species'].str.endswith('i')]
        gk = data[data['Species'].str.endswith('k')]
        gl = data[data['Species'].str.endswith('l')]
        gm = data[data['Species'].str.endswith('m')]
        gn = data[data['Species'].str.endswith('n')]
        gx = data[data['Species'].str.endswith('x')]
        go = data[data['Species'].str.endswith('o')]
        gp = data[data['Species'].str.endswith('p')]
        gr = data[data['Species'].str.endswith('r')]
        gs = data[data['Species'].str.endswith('s')]
        gt = data[data['Species'].str.endswith('t')]
        gu = data[data['Species'].str.endswith('u')]
        gf = data[data['Species'].str.endswith('f')]
        gc = data[data['Species'].str.endswith('c')]
        gy = data[data['Species'].str.endswith('y')]
        gw = data[data['Species'].str.endswith('w')]

        final_data = pd.concat([ga, gb, gg, gd, ge, gz,
                                gh, gq, gi, gk, gl, gm,
                                gn, gx, go, gp, gr, gs,
                                gt, gu, gf, gc, gy, gw],
                                axis=0).reset_index(drop=True)
        
        return final_data