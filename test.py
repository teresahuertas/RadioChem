import pandas as pd
import matplotlib.pyplot as plt
import catalogues as cat
import spectrumfit as sf

#rrls = cat.Catalogue('rrls')
#print(rrls.catalogue)
#ic418 = cat.Catalogue('IC418')
#print(ic418.rrls)
#print(ic418.molecules)
#print(ic418.uf)

files = ['IC418_3mm_tmb_fit_velocity.txt']#, 'IC418_2mm_tmb_fit_velocity.txt']

spectrum = sf.read_spectrum(1, files)
#print(spectrum)
spectrum = sf.create_spectrum(spectrum, 87317.0, 42.59999847412109)
#print(spectrum.unit)

#lines = sf.line_inspector(spectrum, 0.002, 'IC418', 'emission')
#print(lines)

#sf.plot_spectrum(spectrum, lines)

#freq = sf.line_freq_width(87317.0, 60)
#print(freq)

'''plt.plot(spectrum.spectral_axis, spectrum.flux)
plt.xlabel(spectrum.spectral_axis.unit)
plt.ylabel(spectrum.flux.unit)
plt.show()'''