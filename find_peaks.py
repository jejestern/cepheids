import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from sys import argv, exit
import os
from photutils import centroid_1dg, CircularAperture, aperture_photometry, find_peaks
from datetime import *
from astropy.visualization import simple_norm
from astropy.table import Table

# This part takes the argument and saves the folder 
if not len(argv) == 2:
    print("Wrong number of arguments!")
    print("Usage: python Master.py star")
    print("Exiting...")
    exit()

star = argv[1]
files = os.listdir(star)

if star == 'V0473_Lyr':
	x_1 = 2100
	x_2 = 2450
	y_1 = 1700
	y_2 = 2300
	r = 10.
if star == 'RR_Lyr':
	x_1 = 0
	x_2 = 4090
	y_1 = 0
	y_2 = 4090
	r = 10.
if star == 'FF_Aql':
	x_1 = 1700
	x_2 = 2500
	y_1 = 1700
	y_2 = 2500
	r = 10.

for image_name in files[0:1]:
	image = fits.getdata(star+"/"+image_name, ext=0)
	data = image[x_1:x_2, y_1:y_2]
	tbl = find_peaks(data, 300, box_size=50)
	norm = simple_norm(data, 'power', percent=99.9)
	plt.imshow(data, cmap='Greys_r', origin='lower', norm=norm)
	pe = np.array(tbl['peak_value'])
	pe_x = np.array(tbl['x_peak'])
	pe_y = np.array(tbl['y_peak'])
	peaks = np.array((pe_x, pe_y, pe)).T
	peaks = peaks.tolist()
	peaks = sorted(peaks, key=lambda t: t[2], reverse=True)
	positions = ([peaks[0][0], peaks[0][1]], [peaks[1][0], peaks[1][1]]) 
	apertures = CircularAperture(positions, r=r)
	apertures.plot(color ='#0547f9', lw=1.5)
	plt.xlim(0, data.shape[1]-1)
	plt.ylim(0, data.shape[0]-1)
	#plt.savefig(image_name+'.png') # 'Center_picture/'+
	#plt.clf()
	plt.show()
	phot_table = aperture_photometry(data, apertures)
	apertures = np.array(phot_table['aperture_sum'])
	print(apertures)
	
