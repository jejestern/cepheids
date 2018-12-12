import matplotlib.pyplot as plt
from astropy.io import fits
from sys import argv, exit
import os
from astropy.visualization import simple_norm

# This part takes the argument and saves the folder 
if not len(argv) == 2:
    print("Wrong number of arguments!")
    print("Usage: python Master.py file.fit")
    print("Exiting...")
    exit()

image = argv[1]

image = fits.getdata(image, ext=0)
norm = simple_norm(image, 'linear', percent=99.9)
plt.imshow(image, cmap='Greys_r', origin='lower', norm=norm)
plt.show()
