import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from sys import argv, exit
import os
from photutils import CircularAperture, aperture_photometry, find_peaks, CircularAnnulus
from datetime import *
from astropy.table import Table
from astropy.stats import LombScargle 

# This part takes the argument and saves the folder 
if not len(argv) == 2:
    print("Wrong number of arguments!")
    print("Usage: python Master.py star")
    print("Exiting...")
    exit()

star = argv[1]
files = os.listdir(star)

if star == 'RR_Lyr':
	x_1 = 0
	x_2 = 4090
	y_1 = 0
	y_2 = 4090
	r = 8.
	mag_ref = 8.63
	ref = 1
	min_freq = 1.0
	max_freq = 2.0
if star == 'V0473_Lyr':
	x_1 = 2100
	x_2 = 2450
	y_1 = 1700
	y_2 = 2300
	r = 10. 
	mag_ref = 9.16
	ref = 1
	min_freq = 0.3
	max_freq = 1.0
	

time_mag = []
reference_time = datetime(2018, 10, 4, 18, 40, 0, 0)

for image_name in files:
	image = fits.getdata(star+"/"+image_name, ext=0)
	data = image[x_1:x_2, y_1:y_2]
	image_head = fits.getheader(star+"/"+image_name, ext=0)
	tbl = find_peaks(data, 300, box_size=50)
	pe = np.array(tbl['peak_value'])
	pe_x = np.array(tbl['x_peak'])
	pe_y = np.array(tbl['y_peak'])
	peaks = np.array((pe_x, pe_y, pe)).T
	peaks = peaks.tolist()
	peaks = sorted(peaks, key=lambda t: t[2], reverse=True)
	positions = ([peaks[0][0], peaks[0][1]], [peaks[ref][0], peaks[ref][1]])
	"""
	x, y = centroid_1dg(image[x_1:x_2, y_1:y_2])
	#fits.writeto(image_name, image[x_1:x_2, y_1:y_2], overwrite=True)"""
	apertures = CircularAperture(positions, r=r)
	annulus = CircularAnnulus(positions, r_in=12., r_out= 15.)
	apers = [apertures, annulus]
	phot_table = aperture_photometry(data, apers, error=np.sqrt(abs(data)))
	bkg_mean = phot_table['aperture_sum_1']/annulus.area()
	bkg_mean_err = 1/annulus.area() * phot_table['aperture_sum_err_1']
	flux = phot_table['aperture_sum_0'] - bkg_mean * apertures.area()
	flux_err = np.sqrt((phot_table['aperture_sum_err_0'])**2 + (apertures.area()*bkg_mean_err)**2)
	#number_of_counts = np.array(phot_table['aperture_sum']) # The number of counts of the star and the reference star in an cicular area around the star
	#ratio_noc = number_of_counts[0]/number_of_counts[1]
	ratio_flux = flux[0]/flux[1]
	ratio_flux_err = np.sqrt((flux_err[0]/flux[1])**2 + (flux[0]*flux_err[1]/flux[1]**2)**2)
	magnitude = -2.5*np.log10(ratio_flux) + mag_ref # Calculate the apparent magnitude
	magnitude_err = 2.5/(np.log(10)*ratio_flux) * ratio_flux_err
	date = image_head['Date-OBS']
	time = datetime(int(date[:4]), int(date[5:7]), int(date[8:10]), int(date[11:13]), int(date[14:16]), int(date[17:19]), int(float(date[19:22])*10**6))
	time_difference = (float((time - reference_time).total_seconds()))/3600.0/24.0 # Time difference between the reference time and the taking of the image
	time_mag.append((time_difference, magnitude, magnitude_err))
	
times_mag = np.array(sorted(time_mag, key=lambda t: t[0])) 
times = []
time_err = 1/(24*3600)
magnitudes = []
magnitudes_err = []
for i in range(len(files)):
	times.append(times_mag[i][0])
	magnitudes.append(times_mag[i][1])
	magnitudes_err.append(times_mag[i][2])
times = np.array(times)
mean_mag_data = np.mean(magnitudes) # Find the mean magnitude
mean_mag_data_err = np.mean(magnitudes_err)
magnitudes = np.array(magnitudes)
magnitudes_err = np.array(magnitudes_err)

### Calculating the period
#frequency, power = LombScargle(times, magnitudes, magnitudes_err).autopower(minimum_frequency=min_freq, maximum_frequency=max_freq)
frequency = np.linspace(min_freq, max_freq, 100000)
power = LombScargle(times, magnitudes, magnitudes_err).power(frequency)
f = plt.figure()
ax = f.add_subplot(111)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.plot(frequency, power, 'blue', label="periodogram showing possible frequencies")
plt.xlabel('frequencies [1/days]')
plt.ylabel('power')
plt.title('Periodogram '+star)
plt.legend()
plt.savefig('Periodogram_'+star+'.pdf')
plt.clf()
x_peak = np.argmax(power)
freq = frequency[x_peak]
x_L_halb = np.where(power[:x_peak] <= power[x_peak]/2.0)[0][-1]
x_R_halb = np.where(power[:x_peak] <= power[x_peak]/2.0)[0][0]
freq_halb = abs((frequency[x_R_halb]-frequency[x_L_halb])/2.0)
freq_err = freq_halb/np.sqrt(power[x_peak])
Period = 1/ freq
Period_err = 1/freq**2 * freq_err
print('The period of the star is P = ', Period, '+/-', Period_err, 'days.')
prob = LombScargle(times, magnitudes, magnitudes_err).false_alarm_probability(freq)
print('The false alarm probability is ', prob)

### Plotting
f = plt.figure()
ax = f.add_subplot(111)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.errorbar(times, magnitudes, magnitudes_err, time_err, 'r.', label="Magnitudes")
plt.xlabel('time in days [d]')
plt.ylabel('magnitude [m]')
plt.title(star)
plt.legend()
plt.savefig('magnitudes_'+star+'.pdf')
plt.clf()

# Plot the data within one period including the fit
time_new = times
for i in range(len(times)):
	while time_new[i] > Period:
		time_new[i] -= Period

t_fit = np.linspace(0, Period, 1000)
y_fit = LombScargle(times, magnitudes, magnitudes_err).model(t_fit, freq)
f = plt.figure()
ax = f.add_subplot(111)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.errorbar(time_new, magnitudes, magnitudes_err, time_err, 'r.', label="Data points all plotted\nin the same period")
plt.plot(t_fit, y_fit, 'k--', label="Fit of the pulsation of the star")
plt.xlabel('time in days [d]')
plt.ylabel('magnitude [m]')
plt.title(star)
plt.legend()
art1 = []
lgd1 = plt.legend(bbox_to_anchor=(1.03, 1.0))
art1.append(lgd1)
plt.savefig('Period_'+star+'.pdf', additional_artists=art1, bbox_inches='tight')
plt.clf()

### Calculation of the Distance
mean_mag = min(y_fit) + (max(y_fit)-min(y_fit))/2
mean_mag_err = abs(mean_mag-mean_mag_data) + mean_mag_data_err
print("The mean apparent magnitude of", star, "is", mean_mag,"+/-", mean_mag_err)
if star == 'RR_Lyr':
	M = 1.094 + 0.232*-1.16
	M_err = 0.091 + 0.02*-1.16
else:
	M = -2.81*np.log10(Period) - 1.43
	M_err = abs(2.81/(Period*np.log(10))*Period_err)
print('The absolut magnitude of the star is M = ', M, '+/-', M_err)
d = 10**((mean_mag-M+5)/5)
d_err = np.sqrt((np.log(10)/5*10**((mean_mag-M+5)/5)*mean_mag_err)**2 + (np.log(10)/5*10**((mean_mag-M+5)/5)*M_err)**2)
print('The distance is ', d, '+/-', d_err, ' pc.')

