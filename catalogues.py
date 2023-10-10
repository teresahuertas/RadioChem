import os
import numpy as np
import pandas as pd


os.path.abspath(os.getcwd())


class Catalogue:
    """
    Class to read, classify and analyse the obseved lines in PNe

    Attributes
    ----------
    catalogue : pandas dataframe
        Dataframe with the catalogue information
    rrls : pandas dataframe
        Dataframe with the radio recombination line catalogue information sorted by
        element and series
    molecules : pandas dataframe
        Dataframe with the molecular line catalogue information sorted by frequency
    uf : pandas dataframe
        Dataframe with the unidentified line catalogue information sorted by frequency

    Methods
    -------
    read_catalogue_file(name)
        Function to read the catalogue file and return a pandas dataframe
    sort_greek(data)
        Function to sort lines by greek letter
    classify_rrls(data)
        Function to classify the rrls by element and series
    get_rrls(data)
        Function to get the rrls from an astronomical source catalogue
    get_molecules(data)
        Function to get the molecules from an astronomical source catalogue
    get_uf(data)
        Function to get the unidentified lines from an astronomical source catalogue

    Parameters
    ----------
    name : str
        Name of the catalogue file to read
    cat_path : str
        Path to the catalogue files
    results_path : str
        Path to the results folder

    Returns
    -------
    catalogue : pandas dataframe
        Dataframe with the catalogue information
    rrls : pandas dataframe
        Dataframe with the radio recombination line catalogue information sorted by
        element and series
    molecules : pandas dataframe
        Dataframe with the molecular line catalogue information sorted by frequency
    uf : pandas dataframe
        Dataframe with the unidentified line catalogue information sorted by frequency

    Raises
    ------
    FileNotFoundError
        If the file is not found
    Exception
        If there is an error reading the file

    Notes
    -----
    The catalogue files must be in the Source_Catalogues folder and must have the
    following format:
             1         2         3         4         5         6         7         8         9         0
    1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    ----------------------------------------------------------------------------------------------------
    S  # Species        Freq[MHz] Err[MHz] Eup[K]  Gup  Aij[s-1]    Upper level -- Lower level    Origin
    ----------------------------------------------------------------------------------------------------
    C  1 H185\ga         1030.251   0.000     0.0    1  0.00e+00            186 -- 185            rrline

    The columns are:
    S. Status: D for detected, T for tentative, ? for doubtful, F for failed and C for calculated
    #. Number of lines displayed (only makes sense when checking lines with summaryV7multi.class)
    Species. Species name
    Freq[MHz]. Frequency of the line in MHz
    Err[MHz]. Error in the frequency of the line in MHz
    Eup[K]. Upper energy level in Kelvin
    Gup. Upper degeneracy
    Aij[s-1]. Einstein coefficient in s^-1
    Upper level -- Lower level. Upper and lower levels of the transition
    Origin. Origin of the line: rrline for radio recombination line, unknow for UFs and jpl and
            cdms for molecular lines
    """

    # Paths to catalogues defined from the current directory
    cat_path = './Source_Catalogues/'
    results_path = 'RadioChem/Results/'

    # Milimeter bands
    mmbands = np.array([['0.9mm', 277.0e3, 375.0e3],
                        ['1mm', 202.0e3, 274.0e3],
                        ['2mm', 125.0e3, 184.0e3],
                        ['3mm', 73.0e3, 117.0e3],
                        ['7mm', 30.0e3, 50.0e3],
                        ['13mm', 16.0e3, 27.0e3],
                        ['25mm', 8.0e3, 12.0e3],
                        ['5cm', 4.0e3, 8.0e3],
                        ['10cm', 2.0e3, 4.0e3],
                        ['20cm', 1.0e3, 2.0e3]])
    
    # Corresponding IEEE band names
    ieee_bands = np.array([['0.9mm', 'mm'],
                           ['1mm', 'mm'],
                           ['2mm', 'mm'],
                           ['3mm', 'F'],
                           ['7mm', 'Q'],
                           ['13mm', 'K'],
                           ['25mm', 'X'],
                           ['5cm', 'C'],
                           ['10cm', 'S'],
                           ['20cm', 'L']])
    
    # Antennas and observing windows
    telescopes = np.array([['IRAM30m', 225.05e3, 232.837e3],
                           ['IRAM30m', 130.85e3, 139.1e3],
                           ['IRAM30m', 81.5e3, 89.76e3],
                           ['Yebes40m', 31.53e3, 50e3]])

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
                self.uf = self.get_uf(self.catalogue)

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
            # Change data type of columns
            #data['Status'] = data['Status'].astype(str)
            #data['Species'] = data['Species'].astype(str)
            data['Freq[MHz]'] = data['Freq[MHz]'].astype(float)
            #data['Upper'] = data['Upper'].astype(int)
            #data['Lower'] = data['Lower'].astype(int)
            #data['Origin'] = data['Origin'].astype(str)
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
        final_data = self.set_band(final_data)

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
        final_data = self.set_band(final_data)
        
        return final_data.reset_index(drop=True)
    

    def get_uf(self, data):
        """
        Function to get the unidentified lines from an astronomical source catalogue

        Parameters
        ----------
        data : pandas dataframe
            Dataframe with the catalogue information

        Returns
        -------
        final_data : pandas dataframe
            Dataframe with the catalogue information sorted by frequency
        """
        final_data = data[data['Species'].str.startswith('U-')]
        final_data = self.set_band(final_data)
        #final_data = self.set_observed_antenna(final_data)
        
        return final_data.reset_index(drop=True)
    

    @classmethod
    def set_band(cls, data):
        """
        Function to set the frequency band ID to each line
        
        Parameters
        ----------
        data : pandas dataframe
            Dataframe with the catalogue information
            
        Returns
        -------
        data : pandas dataframe
            Dataframe with the catalogue information and the frequency band ID
        """
        # Create a new column called 'Band'
        data['Band'] = None
        # For each line, check the frequency and assign the corresponding band ID in 'Band'
        for i in range(len(cls.mmbands)):
            data.loc[(data['Freq[MHz]'] >= float(cls.mmbands[i,1])) 
                     & (data['Freq[MHz]'] <= float(cls.mmbands[i,2])), 'Band'] = cls.mmbands[i,0]
            if cls.mmbands[i,0] in ['7mm', '13mm', '25mm', '5cm', '10cm', '20cm']:
                data['Band'] = cls.ieee_bands[i,0]
        #for index, row in data.iterrows():
        #    for i in range(len(cls.mmbands)):
                # Check the frequency of the line and assign the corresponding band ID
        #        if cls.mmbands[list(cls.mmbands.keys())[i]][0] <= row['Freq[MHz]'] <= cls.mmbands[list(cls.mmbands.keys())[i]][1]:
        #            data.loc[index, 'Band'] = list(cls.mmbands.keys())[i]
                    # If 'Band' is 7mm, 13mm, 25mm, 5cm, 10cm or 20cm, change the name to the
                    # corresponding IEEE band name
        #            if list(cls.mmbands.keys())[i] in ['7mm', '13mm', '25mm', '5cm', '10cm', '20cm']:
        #                data.loc[index, 'Band'] = cls.ieee_bands[list(cls.mmbands.keys())[i]]
        #            break
        #        else:
        #            continue

        return data
    

    @classmethod
    def set_observed_antenna(cls, data):
        """
        Function to set the antenna used to get the data depending on the observed range

        Parameters
        ----------
        data : pandas dataframe
            Dataframe with the catalogue information

        Returns
        -------
        data : pandas dataframe
            Dataframe with the catalogue information and the antenna used
        """
        data['Telescope'] = None
        # For each line, check the frequency and assign the corresponding antenna
        for index, row in data.iterrows():
            for i in range(len(cls.telescopes)):
                # Check the frequency of the line and assign the corresponding antenna
                if cls.telescopes[list(cls.telescopes.keys())[i]][0] <= row['Freq[MHz]'] <= cls.telescopes[list(cls.telescopes.keys())[i]][1]:
                    data.loc[index, 'Telescope'] = list(cls.telescopes.keys())[i]
                    print(data.loc[index, 'Freq[MHz]'], list(cls.telescopes.keys())[i])
                    break
                else:
                    continue

        return data