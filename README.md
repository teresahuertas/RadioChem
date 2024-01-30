# RadioChem

RadioChem is a python package for the analysis of radiochemistry data. It is designed to combine observational data with rotational transition catalogues and theoretical results. It is currently under development.

## Contents

- catalogues.py: Contains the definition of the class Catalogue as well as the functions to read, write and analyse data from the catalogues.
- spectrumfit.py: Contains the functions to find spectral lines and fit the spectra. It is under development.
- synthetics.py: Contains the functions to work with synthetic spectra. It is under development.

## Usage

```python
import radiochem
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## Examples of use

The following examples show how to use RadioChem to analyse the rotational transitions of a molecule.

### Read a catalogue

```python
from radiochem import catalogues as cat

source = cat.Catalogue(path='path/to/catalogue', name='name_of_source')
```

source.rrls is a pandas dataframe with the radio recombination lines appearing in the catalogue.
source.molecules is a pandas dataframe with the molecules appearing in the catalogue.
source.uf is a pandas dataframe with the unidentified features appearing in the catalogue.

### Read observational spectrum data
``` python
from radiochem import spectrumfit as sf

spectrum = sf.create_spectrum(path='path/to/spectrum/file', name='name_of_spectrum')
```

It is assumed that the spectral data is given in velocity and flux. This function reads data from a file and creates a Spectrum object.

### Find spectral lines

```python
from radiochem import spectrumfit as sf

spectrum = sf.create_spectrum(path='path/to/spectrum/file', name='name_of_spectrum')
lines = sf.line_inspector(spectrum, rms='rms', source='source_name', line_type='line_type')
```

This function finds spectral lines in a spectrum. It returns a QTable with the spectral lines found in the spectrum. The rms parameter is the rms of the spectrum. The source parameter is the name of the source. The line_type parameter is the type of line to be found. It can be 'emission' or 'absorption'.

### Plot a spectrum

```python
from radiochem import spectrumfit as sf

spectrum = sf.create_spectrum(path='path/to/spectrum/file', name='name_of_spectrum')
lines = sf.line_inspector(spectrum, rms='rms', source='source_name', line_type='line_type') # Optional
sf.plot_spectrum(spectrum, lines='detected_lines')
```

This function plots a spectrum. The lines parameter is optional. If it is not given, the spectrum is plotted without the spectral lines. If it is given, the spectrum is plotted with the spectral lines.

### Plot spectrum with lines from a catalogue and synthetic data
``` python
from radiochem import catalogues as cat
from radiochem import spectrumfit as sf
from radiochem import synthetics as synth

source = cat.Catalogue(path='path/to/catalogue', name='name_of_source')
spectrum = sf.create_spectrum(path='path/to/spectrum/file', name='name_of_spectrum')
synthetic_lines = synth.read_synthetic_spectra(path='path/to/synthetic/spectra', name='name_of_synthetic_spectra_file')

synth.plot_synthetic_spectrum('source_name', spectrum, synthetic_lines, source.rrls, source.uf, source.molecules, names=True)
```

This function plots a spectrum with the spectral lines from a catalogue and synthetic data. The names parameter is optional. If it is not given, the spectrum is plotted without the names of rrls, molecules and ufs from the catalogue. If it is given, the spectrum is plotted with the names of rrls, molecules and ufs from the catalogue.