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
            self.catalogue = read_catalogue_file(name)
            if name == 'rrls':
                self.catalogue = self.catalogue
            else:
                self.rrls = 
                self.molecules =
                self.uf =

        return
            