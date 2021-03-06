'''This file contains miscellaneous functions for the catalog pipeline, namely findBox, bboxToDS9, and approx.'''

import numpy as np

def find_box(loop):
   '''
   creates a bounding box for a given contour path
   loop = data['contours'][0][0]['arr'] # outermost contour (for testing)
   '''
   xmax = loop[0]['x']
   ymax = loop[0]['y']
   xmin = loop[0]['x']
   ymin = loop[0]['y']
   for i in loop:
      if i['x']>xmax:
         xmax = i['x']
      elif i['x']<xmin:
         xmin = i['x']
      if i['y']>ymax:
         ymax = i['y']
      elif i['y']<ymin:
         ymin = i['y']
   return [xmax, ymax, xmin, ymin]

def bbox_to_ds9(bbox, img_size):
   '''
   finds the coordinates of the bbox in DS9's system and the input values for drawing a box in DS9
   bbox = tree.value['bbox'] # outermost bbox (for testing)
   '''
   xmax = bbox[0]
   ymax = bbox[1]
   xmin = bbox[2]
   ymin = bbox[3]
   temp = img_size+1-ymax
   ymax = img_size+1-ymin
   ymin = temp
   new_bbox = [xmax, ymax, xmin, ymin]
   ds9_box = [ (xmax+xmin)/2., (ymax+ymin)/2., xmax-xmin, ymax-ymin ]
   return [new_bbox, ds9_box]

def approx(a, b, uncertainty=1e-5):
   '''determines if two floats are approximately equal'''
   return np.abs(a-b) < uncertainty
