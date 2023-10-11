import pandas as pd
import catalogues as cat
import spectrumfit as sf

#rrls = cat.Catalogue('rrls')
#print(rrls.catalogue)
#ic418 = cat.Catalogue('IC418')
#print(ic418.rrls)
#print(ic418.molecules)
#print(ic418.uf)

spectrum = sf.read_spectrum('IC418_3mm_tmb_fit_velocity.txt')
print(spectrum)