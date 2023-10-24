import pandas as pd
import matplotlib.pyplot as plt
import catalogues as cat
#from tqdm import tqdm
import spectrumfit as sf
from astropy import units as u

# Set the number of catalogues to work with
ncat = 3

'''for i in tqdm (range (ncat), desc='Reading catalogues...'):
    rrls = cat.Catalogue('rrls')
    ic418 = cat.Catalogue('IC418')
    ngc7027 = cat.Catalogue('NGC_7027')
print('Catalogues read successfully!')'''
rrls = cat.Catalogue('rrls')
ic418 = cat.Catalogue('IC418')
ngc7027 = cat.Catalogue('NGC_7027')

#print(rrls.catalogue)
#ic418 = cat.Catalogue('IC418')
#print(ic418.rrls)
#print(ngc7027.molecules)
#print(ngc7027.uf)

files = ['IC418_3mm_tmb_fit_velocity.txt']#, 'IC418_2mm_tmb_fit_velocity.txt']
C60CNp = sf.read_synthetic_spectra('C60CN+/C60CN+.dat')
C60CNm = sf.read_synthetic_spectra('C60CN-/C60CN-.dat')
C60CN = sf.read_synthetic_spectra('C60CN0/C60CN.dat')
C60Fep = sf.read_synthetic_spectra('C60Fe+/C60Fe+_eta66.dat')
C60Fe = sf.read_synthetic_spectra('C60Fe0/C60Fe_eta66.dat')
C70p = sf.read_synthetic_spectra('C70+/C70+_2.73K.txt')
#print(C60CNp)

spectrum = sf.read_spectrum(1, files)
#print(spectrum)
spectrum = sf.create_spectrum(spectrum, 87317.0, 42.59999847412109, 12.0)
#print(spectrum.unit)

#lines = sf.line_inspector(spectrum, 0.002, 'IC418', 'emission')
# Save lines to a file
#with open('lines.txt', 'w') as f:
#    for line in lines:
#        f.write(f'{line}\n')
#print(lines)

#sf.plot_spectrum(spectrum, lines)
#sf.plot_synthetic_spectrum(spectrum, C60CNp)

#freq = sf.line_freq_width(87317.0, 60)
#print(freq)

'''plt.plot(spectrum.spectral_axis, spectrum.flux)
plt.xlabel(spectrum.spectral_axis.unit)
plt.ylabel(spectrum.flux.unit)'''
# Plot synthetic spectrum
fig, ax = plt.subplots()
ax.plot(spectrum.spectral_axis, spectrum.flux)
ax.hlines(0.0, min(spectrum.spectral_axis/u.MHz), 
              max(spectrum.spectral_axis/u.MHz), color='black')
ax.vlines(ic418.rrls['Freq[MHz]'], 0, max(spectrum.flux/u.K), 
          colors='b', label='rrls')
ax.vlines(ic418.uf['Freq[MHz]'], 0, max(spectrum.flux/u.K), 
          colors='g', label='uf')
ax.vlines(C60CNp['Freq[MHz]'], 0, max(spectrum.flux/u.K), 
          colors='r', label='C60CN+')
ax.vlines(C60CNm['Freq[MHz]'], 0, max(spectrum.flux/u.K),
            colors='c', label='C60CN-')
ax.vlines(C60CN['Freq[MHz]'], 0, max(spectrum.flux/u.K),
            colors='m', label='C60CN')
ax.vlines(C60Fep['Freq[MHz]'], 0, max(spectrum.flux/u.K),
            colors='brown', label='C60Fe+')
ax.vlines(C60Fe['Freq[MHz]'], 0, max(spectrum.flux/u.K),
            colors='y', label='C60Fe')
ax.vlines(C70p['Freq[MHz]'], 0, max(spectrum.flux/u.K),
            colors='k', label='C70+')
# Add legend outside the plot
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
# Set x-axis limits
ax.set_xlim(min(spectrum.spectral_axis/u.MHz), 
            max(spectrum.spectral_axis/u.MHz))
ax.set_xlabel(spectrum.spectral_axis.unit)
ax.set_ylabel(spectrum.flux.unit)
plt.show()