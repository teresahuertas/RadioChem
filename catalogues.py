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
                self.catalogue = self.classify_rrls(self.catalogue)
            else:
                self.rrls = self.get_rrls(self.catalogue)
                self.molecules = self.get_molecules(self.catalogue)
                #self.uf =

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
    

    def classify_rrls(self, data):
        """
        Function to classify the rrls by element and series

        Parameters
        ----------
        data : pandas dataframe
            Dataframe with the catalogue information

        Returns
        -------
        final_data : pandas dataframe
            Dataframe with the catalogue information sorted by element and series
        """
        # Get hydrogen lines
        # They start by H and not by He or 3He
        Hydrogen = data[data['Species'].str.startswith('H') 
                        & ~data['Species'].str.startswith('He') 
                        & ~data['Species'].str.startswith('3He')]
        Deuterium = data[data['Species'].str.startswith('D')]
        TriHelium = data[data['Species'].str.startswith('3He') 
                         & ~data['Species'].str.startswith('3HeII')]
        Helium = data[data['Species'].str.startswith('He')]
        Carbon = data[data['Species'].str.startswith('C') 
                      & ~data['Species'].str.startswith('CII') 
                      & ~data['Species'].str.startswith('CIII')]
        TriHeliumII = data[data['Species'].str.startswith('3HeII')]
        HeliumII = data[data['Species'].str.startswith('HeII')]
        CarbonII = data[data['Species'].str.startswith('CII') 
                        & ~data['Species'].str.startswith('CIII')]
        CarbonIII = data[data['Species'].str.startswith('CIII')]
        OxygenIII = data[data['Species'].str.startswith('OIII')]

        # Sort lines by greek letter
        Hydrogen = self.sort_greek(Hydrogen)
        Deuterium = self.sort_greek(Deuterium)
        TriHelium = self.sort_greek(TriHelium)
        Helium = self.sort_greek(Helium)
        Carbon = self.sort_greek(Carbon)
        TriHeliumII = self.sort_greek(TriHeliumII)
        HeliumII = self.sort_greek(HeliumII)
        CarbonII = self.sort_greek(CarbonII)
        CarbonIII = self.sort_greek(CarbonIII)
        OxygenIII = self.sort_greek(OxygenIII)

        final_data = pd.concat([Hydrogen, Deuterium, TriHelium, Helium, 
                                Carbon, TriHeliumII, HeliumII, CarbonII, 
                                CarbonIII, OxygenIII], axis=0).reset_index(drop=True)

        return final_data
    

    def get_rrls(self, data):
        """
        Function to get the rrls from an astronomical source catalogue

        Parameters
        ----------
        data : pandas dataframe
            Dataframe with the catalogue information

        Returns
        -------
        self.rrls : pandas dataframe
            Dataframe with the catalogue information sorted by element and series
        """
        final_data = data[data['Origin'] == 'rrline']
        final_data = self.classify_rrls(final_data)

        return final_data.reset_index(drop=True)
    

    def get_molecules(self, data):
        """
        Function to get the molecules from an astronomical source catalogue

        Parameters
        ----------
        data : pandas dataframe
            Dataframe with the catalogue information

        Returns
        -------
        self.molecules : pandas dataframe
            Dataframe with the catalogue information sorted by element and series
        """
        final_data = data[~data['Origin'].isin(['rrline', 'unknow'])]
        
        return final_data.reset_index(drop=True)