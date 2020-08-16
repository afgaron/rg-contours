#custom modules for the RGZ catalog
import processing as p #contains functions that process the data

# Example inputs
fits_file = '/home/garo0040/Documents/rg-contours/FIRSTJ124646.7+224554.fits'
contour_file = '/home/garo0040/Documents/rg-contours/52af7d53eb9a9b05ef000809.json'
host_ra = 191.6946556
host_dec = 22.7615058

source = {'host_data':{'ra':host_ra, 'dec':host_dec},
		  'fits_file':fits_file,
		  'contour_file':contour_file}

# Read radio contour and flux data from file
radio_data = p.get_radio(source)
source['radio'] = radio_data

# Calculate bending and orientation from contours
bending_data = p.get_bending(source)
source['bending'] = bending_data