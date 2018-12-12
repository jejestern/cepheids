from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
import os
from astroscrappy import detect_cosmics
# This package is from https://github.com/astropy/astroscrappy
from astropy.stats import sigma_clipped_stats, SigmaClip
from photutils import make_source_mask, Background2D, MedianBackground


# This part takes the argument and saves the folder 
if not len(argv) == 2:
    print("Wrong number of arguments!")
    print("Usage: python Master.py base_directory_name")
    print("Exiting...")
    exit()

input_filename = argv[1]

star_names = ['Gam_Lyr', 'RR_Lyr', 'V0473_Lyr_1', 'V0473_Lyr_05', 'Zeta_Aql'] # 'FF_Aql_2', 'FF_Aql_02'

# Import the flats and the darks of the flats and c
Dark_for_Flat = os.listdir(input_filename+"/Dark_for_Flat")
dark_flat = []
for item in Dark_for_Flat:
	dark_flat.append(fits.getdata(input_filename+"/Dark_for_Flat/"+item, ext=0))
dark_flat = np.median(dark_flat, axis=0)

Flat = os.listdir(input_filename+"/Flat")
flat = []
for item in Flat:
	flat.append(fits.getdata(input_filename+"/Flat/"+item, ext=0) - dark_flat)

# Creation of the masterflat and its mean
masterflat = np.median(flat, axis=0)
masterflat_mean = np.mean(masterflat)

# We delete the information on the flat and dark_flat
del Dark_for_Flat
del Flat 
del flat
del dark_flat

for star in star_names:
	print(star)
	Dark_Science = os.listdir(input_filename+"/"+star+"_Dark")
	dark_science = []
	for item in Dark_Science:
		dark_science.append(fits.getdata(input_filename+"/"+star+"_Dark/"+item, ext=0))
	dark_science = np.median(dark_science, axis=0) # Masterdark of science
	del Dark_Science
	Science = os.listdir(input_filename+"/"+star)
	for Sci in Science:
		print(Sci)
		science = fits.getdata(input_filename+"/"+star+"/"+Sci, ext=0)
		# Calculations of the final science
		science = (np.array(science) - dark_science)/(masterflat/masterflat_mean)
		head = fits.getheader(input_filename+"/"+star+"/"+Sci, ext=0)
		mask, science = detect_cosmics(indat=science, inmask=None, sigclip=4.0, sigfrac=0.3, objlim=2000, pssl=0.0, gain=head['EGAIN'], readnoise=9.0, satlevel=100000.0, niter=5, sepmed=True, cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=False)
		# Camera specific values http://diffractionlimited.com/wp-content/uploads/2016/08/stx16803_specs_7_12_11.pdf
		
		mask = (science == 0)
		bkg = Background2D(science, (250,250), filter_size=(10,10), sigma_clip = SigmaClip(sigma=3., iters=10), bkg_estimator = MedianBackground(), mask=mask)
				
		fits.writeto(star+'/'+Sci, science - bkg.background* ~mask, head, overwrite=True)
		del mask
		del bkg
	
	

	

