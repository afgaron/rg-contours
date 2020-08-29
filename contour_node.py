'''
This file contains an implementation of a contour tree object, designed for measuring the total and peak fluxes of arbitrary radio galaxies. Each Node contains a contour and links to the immediately interior contours.
The total radio flux through a contour, as well as radio peaks, can be calculated from the FITS file plus the contour pixels.
Disclaimer: it is assumed that the FITS image is square, and the loop in get_total_flux() will fail if it's not.
'''

import numpy as np
from astropy.io import fits
from astropy import wcs
from scipy.special import erfinv
from matplotlib import path

# Custom modules
import misc_functions as fn # Contains miscellaneous helper functions

class Node(object):
    '''Tree implementation for contours'''
    
    def __init__(self, value=None, contour=None, fits_loc=None, img=None, header=None, w=None, sigma_mJy_beam=0):
        '''Tree initializer'''
        self.value = value # Contour curve and level data
        self.children = [] # Next contour curves contained within this one
        if fits_loc is not None:
            self.get_fits(fits_loc)
        else:
            self.img = img # FITS data as an array
            self.img_size = int(img.shape[0]) # Size in pixels of FITS data
            self.header = header
            self.w = w # WCS converter object
        self.beam_area_arcsec2 = 1.44/4*np.pi*self.header['BMAJ']*self.header['BMIN']*3600*3600 # FWHM ellipse
        self.pixel_area_arcsec2 = wcs.utils.proj_plane_pixel_area(self.w)*3600*3600 # Arcsecond^2
        if contour is not None:
            mad2sigma = np.sqrt(2)*erfinv(2*0.75-1) # Conversion factor
            self.sigma_mJy_beam = 1000*(contour[0]['level']/3) / mad2sigma # Standard deviation of flux density measurements
            for i in contour:
                self.insert(Node(value=i, img=self.img, header=self.header, w=self.w, sigma_mJy_beam=self.sigma_mJy_beam))
            vertices = []
            for pos in contour[0]['arr']:
                vertices.append([pos['x'], pos['y']])
            self.path_outline = path.Path(vertices) # self.path_outline is a Path object tracing the contour
            self.get_total_flux() # self.flux_mJy and self.flux_err_mJy are the total integrated flux and error, respectively; also sets self.area_arcsec2, in arcsec^2
            self.get_peaks() # self.peaks is list of dicts of peak fluxes and locations
        else:
            self.sigma_mJy_beam = sigma_mJy_beam
            self.path_outline = None
            self.flux_mJy = 0
            self.flux_err_mJy = 0
            self.peaks = []
    
    def insert(self, new_node):
        '''Insert a contour node'''
        if self.value is None: # Initialize the root with the outermost contour
            self.value = new_node.value
        elif self.value == new_node.value: # No duplicate contours
            return
        else:
            if new_node.value['k'] == self.value['k'] + 1: # Add a contour one level higher as a child
                self.children.append(new_node)
            elif new_node.value['k'] <= self.value['k']: # If a contour of lower level appears, something went wrong
                raise RuntimeError('Inside-out contour')
            else: # Otherwise, find the next level that has a bounding box enclosing the new contour
                inner = fn.find_box(new_node.value['arr'])
                for child in self.children:
                    outer = fn.find_box(child.value['arr'])
                    if outer[0]>inner[0] and outer[1]>inner[1] and outer[2]<inner[2] and outer[3]<inner[3]:
                        child.insert(new_node)
    
    def check(self):
        '''Manually check the topology of the tree by printing level numbers and bounding boxes to screen (for testing only)'''
        if self.value is None:
            print 'Empty'
        else:
            print 'Level {}: {}'.format(self.value['k'], fn.find_box(self.value['arr']))
            if self.children == []:
                print 'End'
            else:
                for child in self.children:
                    child.check()
    
    def get_fits(self, fits_loc):
        '''Read FITS data from file'''
        self.img = fits.getdata(fits_loc, 0) # Imports data as array
        self.img[np.isnan(self.img)] = 0 # Sets NANs to 0
        self.img_size = int(self.img.shape[0])
        self.header = fits.getheader(fits_loc, 0)
        self.w = wcs.WCS(self.header) # Gets pixel-to-WCS conversion from header
        return self.img
    
    def get_total_flux(self):
        '''Find the total integrated flux of the component and its error'''
        flux_density_Jy_beam = 0
        n = 0
        for i in range(self.img_size-1):
            for j in range(self.img_size-1):
                if self.contains([i+1, self.img_size-j]):
                    flux_density_Jy_beam += self.img[j][i]
                    n += 1
        self.area_arcsec2 = n*self.pixel_area_arcsec2
        flux_density_err_mJy_beam = np.sqrt(n)*self.sigma_mJy_beam
        self.flux_mJy = flux_density_Jy_beam*1000*self.pixel_area_arcsec2/self.beam_area_arcsec2
        self.flux_err_mJy = flux_density_err_mJy_beam*self.pixel_area_arcsec2/self.beam_area_arcsec2
        return [self.flux_mJy, self.flux_err_mJy]
    
    def get_peaks(self, peak_list=None):
        '''Finds the peak values (in mJy) and locations (in DS9 pixels) and return as dict'''
        if peak_list is None:
            peak_list = []
        if self.children == []:
            bbox = fn.bbox_to_ds9(fn.find_box(self.value['arr']), self.img_size)[0] # Bounding box of innermost contour
            flux_density_Jy_beam = self.img[ int(bbox[3]-1):int(bbox[1]+1), int(bbox[2]-1):int(bbox[0]+1) ].max() # Peak flux in bbox, with 1 pixel padding
            loc_pix = np.where(self.img == flux_density_Jy_beam) # Location in pixels
            loc_wcs = self.w.wcs_pix2world( np.array( [[loc_pix[1][0]+1, loc_pix[0][0]+1]] ), 1) # Location in RA and Dec
            peak = dict(ra=loc_wcs[0][0], dec=loc_wcs[0][1], flux=flux_density_Jy_beam*1000)
            peak_list.append(peak)
        else:
            for child in self.children:
                child.get_peaks(peak_list)
        self.peaks = peak_list
        return self.peaks
    
    def contains(self, point):
        '''Returns 1 if point is within the contour, returns 0 otherwise or if there is no contour data'''
        if self.path_outline is not None:
            return self.path_outline.contains_point(point)
        else:
            return 0
