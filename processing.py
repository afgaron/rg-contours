from astropy import coordinates as coord, units as u
import numpy as np

#custom modules for the RGZ catalog
import misc_functions as fn #contains miscellaneous helper functions
import contour_node as c #contains Node class

def getRadio(data, fits_loc, consensusObject):
    '''
    calculates all of the radio parameters from the fits file
    data is a JSON object containing the contours
    fits_loc is a path to the FITS file
    '''
    
    #create list of trees, each containing a contour and its contents
    contourTrees = []
    consensusBboxes = consensusObject['bbox']
    for contour in data['contours']:
        for bbox in consensusBboxes:
            if fn.approx(contour[0]['bbox'][0], bbox[0]) and fn.approx(contour[0]['bbox'][1], bbox[1]) and \
               fn.approx(contour[0]['bbox'][2], bbox[2]) and fn.approx(contour[0]['bbox'][3], bbox[3]):
                tree = c.Node(contour=contour, fits_loc=fits_loc)
                contourTrees.append(tree)
    
    #get component fluxes and sizes
    components = []
    for tree in contourTrees:
        bboxP = fn.bboxToDS9(fn.findBox(tree.value['arr']), tree.imgSize)[0] #bbox in DS9 coordinate pixels
        bboxCornersRD = tree.w.wcs_pix2world( np.array( [[bboxP[0],bboxP[1]], [bboxP[2],bboxP[3]] ]), 1) #two opposite corners of bbox in ra and dec
        raRange = [ min(bboxCornersRD[0][0], bboxCornersRD[1][0]), max(bboxCornersRD[0][0], bboxCornersRD[1][0]) ]
        decRange = [ min(bboxCornersRD[0][1], bboxCornersRD[1][1]), max(bboxCornersRD[0][1], bboxCornersRD[1][1]) ]
        pos1 = coord.SkyCoord(raRange[0], decRange[0], unit=(u.deg, u.deg))
        pos2 = coord.SkyCoord(raRange[1], decRange[1], unit=(u.deg, u.deg))
        extentArcsec = pos1.separation(pos2).arcsecond
        solidAngleArcsec2 = tree.areaArcsec2
        components.append({'flux':tree.fluxmJy, 'flux_err':tree.fluxErrmJy, 'angular_extent':extentArcsec, 'solid_angle':solidAngleArcsec2, \
                           'ra_range':raRange, 'dec_range':decRange})
    
    #adds up total flux of all components
    totalFluxmJy = 0
    totalFluxErrmJy2 = 0
    for component in components:
        totalFluxmJy += component['flux']
        totalFluxErrmJy2 += np.square(component['flux_err'])
    totalFluxErrmJy = np.sqrt(totalFluxErrmJy2)
    
    #finds total area enclosed by contours in arcseconds
    totalSolidAngleArcsec2 = 0
    for component in components:
        totalSolidAngleArcsec2 += component['solid_angle']
    
    #find maximum extent of component bboxes in arcseconds
    maxAngularExtentArcsec = 0
    if len(components)==1:
        maxAngularExtentArcsec = components[0]['angular_extent']
    else:
        for i in range(len(components)-1):
            for j in range(1,len(components)-i):
                corners1 = np.array([ [components[i]['ra_range'][0], components[i]['dec_range'][0]], \
                                      [components[i]['ra_range'][0], components[i]['dec_range'][1]], \
                                      [components[i]['ra_range'][1], components[i]['dec_range'][0]], \
                                      [components[i]['ra_range'][1], components[i]['dec_range'][1]] ])
                corners2 = np.array([ [components[i+j]['ra_range'][0], components[i+j]['dec_range'][0]], \
                                      [components[i+j]['ra_range'][0], components[i+j]['dec_range'][1]], \
                                      [components[i+j]['ra_range'][1], components[i+j]['dec_range'][0]], \
                                      [components[i+j]['ra_range'][1], components[i+j]['dec_range'][1]] ])
                pos1 = coord.SkyCoord(corners1.T[0], corners1.T[1], unit=(u.deg, u.deg))
                pos2 = coord.SkyCoord(corners2.T[0], corners2.T[1], unit=(u.deg, u.deg))
                angularExtentArcsec = pos1.separation(pos2).arcsecond
                maxAngularExtentArcsec = max(np.append(angularExtentArcsec, maxAngularExtentArcsec))
    
    #combine all peaks into single list
    peakList = []
    for tree in contourTrees:
        for peak in tree.peaks:
            peakList.append(peak)
    peakFluxErrmJy = contourTrees[0].sigmamJy

    #find center of radio source
    raMin, raMax, decMin, decMax = np.inf, 0, np.inf, 0
    for comp in components:
        if comp['ra_range'][0] < raMin:
            raMin = comp['ra_range'][0]
        if comp['ra_range'][1] > raMax:
            raMax = comp['ra_range'][1]
        if comp['dec_range'][0] < decMin:
            decMin = comp['dec_range'][0]
        if comp['dec_range'][1] > decMax:
            decMax = comp['dec_range'][1]
    meanRa = (raMax+raMin)/2.
    meanDec = (decMax+decMin)/2.
    
    radio_data = {'radio':{'total_flux':totalFluxmJy, 'total_flux_err':totalFluxErrmJy, 'outermost_level':data['contours'][0][0]['level']*1000, \
                           'number_components':len(contourTrees), 'number_peaks':len(peakList), 'max_angular_extent':maxAngularExtentArcsec, \
                           'total_solid_angle':totalSolidAngleArcsec2, 'peak_flux_err':peakFluxErrmJy, 'peaks':peakList, 'components':components, \
                           'ra':meanRa, 'dec':meanDec}}
    
    return radio_data
