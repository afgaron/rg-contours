# Custom modules
import processing as p # Contains functions that process the data

# Example inputs
fits_file = 'FIRSTJ124646.7+224554.fits'
contour_file = '52af7d53eb9a9b05ef000809.json'
host_ra = 191.6946556
host_dec = 22.7615058

source = {'host_data':{'ra':host_ra, 'dec':host_dec},
		  'fits_file':fits_file,
		  'contour_file':contour_file}

# Measure the radio peaks, integrated flux, and angular size of the radio galaxy
# Required inputs: fits_file and contour_file
radio_data = p.get_radio(source)
source['radio'] = radio_data

# Calculate bending and orientation from contours
# Required inputs: host_data, fits_file, contour_file, and source['radio']['number_peaks'] (generated by p.get_radio())
bending_data = p.get_bending(source)
source['bending'] = bending_data
