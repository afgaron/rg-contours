from astropy import coordinates as coord, units as u
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import interp1d

# Custom modules
import misc_functions as fn # Contains miscellaneous helper functions
import contour_node as c # Node class for processing radio fluxes
import contour_path_object as cpo # Node class for processing bending information

def get_radio(source):
	'''calculates all of the radio parameters from the fits file'''
	
	# Read radio contour from file
	with open(source['contour_file'], 'r') as cf:
		data = json.load(cf)
	
	# Create list of trees, each containing a contour and its contents
	contourTrees = []
	consensusBboxes = [contour[0]['bbox'] for contour in data['contours']]
	for contour in data['contours']:
		for bbox in consensusBboxes:
			if fn.approx(contour[0]['bbox'][0], bbox[0]) and fn.approx(contour[0]['bbox'][1], bbox[1]) and \
			   fn.approx(contour[0]['bbox'][2], bbox[2]) and fn.approx(contour[0]['bbox'][3], bbox[3]):
				tree = c.Node(contour=contour, fits_loc=source['fits_file'])
				contourTrees.append(tree)
	
	# Get component fluxes and sizes
	components = []
	for tree in contourTrees:
		bboxP = fn.bboxToDS9(fn.findBox(tree.value['arr']), tree.imgSize)[0] # Bounding box in DS9 coordinate pixels
		bboxCornersRD = tree.w.wcs_pix2world( np.array( [[bboxP[0],bboxP[1]], [bboxP[2],bboxP[3]] ]), 1) # Two opposite corners of bbox in RA and Dec
		raRange = [ min(bboxCornersRD[0][0], bboxCornersRD[1][0]), max(bboxCornersRD[0][0], bboxCornersRD[1][0]) ]
		decRange = [ min(bboxCornersRD[0][1], bboxCornersRD[1][1]), max(bboxCornersRD[0][1], bboxCornersRD[1][1]) ]
		pos1 = coord.SkyCoord(raRange[0], decRange[0], unit=(u.deg, u.deg))
		pos2 = coord.SkyCoord(raRange[1], decRange[1], unit=(u.deg, u.deg))
		extentArcsec = pos1.separation(pos2).arcsecond
		solidAngleArcsec2 = tree.areaArcsec2
		components.append({'flux':tree.fluxmJy, 'flux_err':tree.fluxErrmJy, 'angular_extent':extentArcsec, 'solid_angle':solidAngleArcsec2, \
						   'ra_range':raRange, 'dec_range':decRange})
	
	# Adds up total flux of all components
	totalFluxmJy = 0
	totalFluxErrmJy2 = 0
	for component in components:
		totalFluxmJy += component['flux']
		totalFluxErrmJy2 += np.square(component['flux_err'])
	totalFluxErrmJy = np.sqrt(totalFluxErrmJy2)
	
	# Finds total area enclosed by contours in arcseconds
	totalSolidAngleArcsec2 = 0
	for component in components:
		totalSolidAngleArcsec2 += component['solid_angle']
	
	# Find maximum extent of component bboxes in arcseconds
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
	
	# Combine all peaks into single list
	peakList = []
	for tree in contourTrees:
		for peak in tree.peaks:
			peakList.append(peak)
	peakFluxErrmJy = contourTrees[0].sigmamJy
	
	# Find center of radio source
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
	
	radio_data = {'total_flux':totalFluxmJy, 'total_flux_err':totalFluxErrmJy, 'outermost_level':data['contours'][0][0]['level']*1000, \
				  'number_components':len(contourTrees), 'number_peaks':len(peakList), 'max_angular_extent':maxAngularExtentArcsec, \
				  'total_solid_angle':totalSolidAngleArcsec2, 'peak_flux_err':peakFluxErrmJy, 'peaks':peakList, 'components':components, \
				  'ra':meanRa, 'dec':meanDec}}
	
	return radio_data

def get_bending(source):
	'''
	Calculate all the bending parameters that don't depend on the cluster
	'''
	
	peak_count = source['radio']['number_peaks']
	assert (peak_count in [2, 3]), 'Not a valid morphology'
	
	# Get pixel-to-WCS conversion
	w = wcs.WCS(fits.getheader(source['fits_file'], 0))
	
	# Get the location of the source
	ir = coord.SkyCoord(source['host_data']['ra'], source['host_data']['dec'], unit=(u.deg,u.deg), frame='icrs')
	
	ir_pos = w.wcs_world2pix(np.array([[ir.ra.deg,ir.dec.deg]]), 1)
	peaks = source['radio']['peaks']
	peak_pos = w.wcs_world2pix(np.array([ [peak['ra'],peak['dec']] for peak in peaks ]), 1)
	
	# Get image parameters for this source
	with open(source['contour_file'], 'r') as cf:
		data = json.load(cf)
	
	contour_tree = get_contours(w, ir_pos, peak_pos, data, peak_count)
	peaks = get_global_peaks(w, peak_pos, peaks, contour_tree)
	if len(peaks) != 2:
		output("%s didn't have 2 tails" % source['zooniverse_id'])
		return
	
	contour_list = [child.path for child in contour_tree.children if any(child.contains(peak_pos))]
	
	# Using the 'contour' method
	bending_angles = get_angles(w, ir, 'contour', contour_list)
	tail_lengths = get_tail_lengths(w, ir, 'contour', contour_list)
	ratios = peak_edge_ratio(w, ir, peaks, tail_lengths_apparent)
	asymmetry = ratios[1]/ratios[0]
	using_contour = {'tail_deg_0':tail_lengths[0], 'tail_deg_1':tail_lengths[1], 'size_deg':sum(tail_lengths), 'ratio_0':ratios[0], 'ratio_1':ratios[1], 'asymmetry':max(asymmetry,1./asymmetry)}
	using_contour.update(bending_angles)
	for key in using_contour.keys():
		if type(using_contour[key]) is coord.angles.Angle:
			using_contour[key] = using_contour[key].deg
	
	# Using the 'peak' method
	bending_angles = get_angles(w, ir, 'peak', peaks)
	tail_lengths = get_tail_lengths(w, ir, 'peak', contour_list, peaks)
	ratios = peak_edge_ratio(w, ir, peaks, tail_lengths)
	asymmetry = ratios[1]/ratios[0]
	using_peaks = {'tail_deg_0':tail_lengths[0], 'tail_deg_1':tail_lengths[1], 'size_deg':sum(tail_lengths), 'ratio_0':ratios[0], 'ratio_1':ratios[1], 'asymmetry':max(asymmetry,1./asymmetry)}
	using_peaks.update(bending_angles)
	for key in using_peaks.keys():
		if type(using_peaks[key]) is coord.angles.Angle:
			using_peaks[key] = using_peaks[key].deg
	
	morphology = 'double' if peak_count == 2 else 'triple'
	entry = {'morphology':morphology, 'using_contour':using_contour, 'using_peaks':using_peaks}
	
	return entry

def get_contours(w, ir_pos, peak_pos, data, peak_count):
	'''
	Returns a list of Path objects corresponding to each outer contour in the data, in RA and dec coordinates
	Removes outer layers until there are two components and a disjoint IR
	Removes any contours that aren't in this source
	'''
	assert (peak_count in [2,3]), 'Not a valid morphology'
	
	# Assemble the contour trees
	contour_trees = []
	for contour in data['contours']:
		tree = cpo.Node(w, contour=contour)
		contour_trees.append(tree)
	
	# Remove each contour that doesn't contain a peak from this source
	contains_peak = []
	for ix, tree in enumerate(contour_trees):
		if any(tree.contains(peak_pos)):
			contains_peak.append(ix)
	
	contour_trees[:] = [tree for ix, tree in enumerate(contour_trees) if ix in contains_peak]
	
	# Combine the entire source into a single tree
	value_at_inf = {'arr':[{'x':-1,'y':-1}, {'x':w._naxis1+1,'y':-1}, {'x':w._naxis1+1,'y':w._naxis2+1}, {'x':-1,'y':w._naxis2+1}, {'x':-1,'y':-1}], 'k':-1}
	source_tree = cpo.Node(w)
	source_tree.insert(cpo.Node(w, value=value_at_inf))
	source_tree.children = contour_trees
	
	# Remove the BCG source if it's a triple
	if peak_count == 3:
		source_tree.remove_triple_center(ir_pos, peak_pos)
	
	# Increase the contour level until the IR position is outside all the contours
	roots = []
	source_tree.get_equal_disjoint(ir_pos, roots)
	source_tree.children = roots
	
	return source_tree

def get_global_peaks(w, peak_pos, peaks, contour_tree):
	'''
	Determines the position of the global maximum for each component in the contour
	'''
	global_peaks = []
	for child in contour_tree.children:
		global_peak = {'flux':0}
		for peak in [peaks[ix] for ix, elem in enumerate(child.contains(peak_pos)) if elem]:
			if peak['flux'] > global_peak['flux']:
				global_peak = peak
		if global_peak['flux'] > 0:
			global_peaks.append(global_peak)
	return global_peaks

def get_angles(w, ir, method, method_data):
	'''
	Determines the opening angle of the radio tail and the position angle of the angle bisector
	Method:
		Contour: from the IR position to the most distant part of the radio contour (data is contour_list)
		Peak: from the IR position to the peak of each component (data is source['radio']['peaks'])
	'''
	assert (method in ['contour', 'peak']), 'Not a valid method'
	pos_angle_0 = get_pos_angle(w, ir, method, method_data[0])
	pos_angle_1 = get_pos_angle(w, ir, method, method_data[1])
	opening_angle = np.abs(pos_angle_1-pos_angle_0).wrap_at(2*np.pi*u.rad)
	bending_angle = coord.Angle(np.abs(np.pi*u.rad - opening_angle))
	bisector = (pos_angle_1+pos_angle_0)/2.
	if np.abs(bisector-pos_angle_0) > np.pi/2*u.rad:
		bisector += np.pi*u.rad
	bending_angles = {'pos_angle_0':pos_angle_0, 'pos_angle_1':pos_angle_1, 'bending_angle':bending_angle, 'bisector':bisector.wrap_at(2*np.pi*u.rad)}
	return bending_angles

def get_pos_angle(w, ir, method, method_data):
	'''
	Determines the position angle between the IR position and the given comparison object
	Method:
		Contour: from the IR position to the most distant part of the radio contour (data is contour_list)
		Peak: from the IR position to the peak of each component (data is source['radio']['peaks'])
	'''
	if method == 'contour':
		contour_sky = coord.SkyCoord(w.wcs_pix2world(method_data.vertices,1), unit=(u.deg,u.deg), frame='icrs')
		separation = ir.separation(contour_sky)
		pos_angle = ir.position_angle(contour_sky)[separation==np.max(separation)][0]
	elif method == 'peak':
		pos_angle = ir.position_angle(coord.SkyCoord(method_data['ra'], method_data['dec'], unit=(u.deg,u.deg), frame='icrs'))
	return pos_angle

def get_tail_lengths(w, ir, method, contour_list, peaks=None):
	'''
	Determines angular separation between the IR position and the given comparison object
	Method:
		Contour: from the IR position to the most distant part of the radio contour
		Peak: from the IR position to the peak of the component
	'''
	tail_lengths = []
	if method == 'contour':
		for contour in contour_list:
			contour_sky = coord.SkyCoord(w.wcs_pix2world(contour.vertices,1), unit=(u.deg,u.deg), frame='icrs')
			separation = ir.separation(contour_sky)
			tail_lengths.append(np.max(separation))
	elif method == 'peak':
		assert (peaks is not None), 'No radio peaks provided'
		for contour, peak in zip(contour_list, peaks):
			tail_lengths.append(get_colinear_separation(w, ir, peak, contour))
	return tail_lengths

def get_colinear_separation(w, ir, peak, contour):
	'''
	Finds the distance from the host to the edge of the contour, passing through the peak
	'''
	
	ir_pos = w.wcs_world2pix(np.array([[ir.ra.deg,ir.dec.deg]]), 1)[0]
	peak_pos = w.wcs_world2pix(np.array([[peak['ra'], peak['dec']]]), 1)[0]
	
	# Extrapolate the line connecting the peak to the IR position
	slope = (peak_pos[1]-ir_pos[1])/(peak_pos[0]-ir_pos[0])
	extrap_pos = ir_pos + w._naxis1*np.array([1.,slope])
	extrap_neg = ir_pos - w._naxis1*np.array([1.,slope])
	
	# Split the contours into well-behaved functions
	# Roll the array until the first index is the minimum value
	x, y = contour.vertices.T
	xmin_loc = np.where(x==min(x))[0][0]
	x_rot = np.append(np.roll(x[:-1], len(x)-xmin_loc-1), min(x))
	y_rot = np.append(np.roll(y[:-1], len(x)-xmin_loc-1), y[xmin_loc])
	
	# Find where the contour doubles back on itself along the x-axis
	m_sign = np.sign(x_rot[1:]-x_rot[:-1])
	roots = np.where(m_sign[1:] - m_sign[:-1] != 0)[0] + 1
	limits = np.concatenate(([0], roots, [len(x_rot)-1]))
	
	# Split the contours at the double-back positions
	domains = []
	ranges = []
	for ix in range(len(limits)-1):
		domains.append(x_rot[limits[ix]:limits[ix+1]+1])
		ranges.append(y_rot[limits[ix]:limits[ix+1]+1])
	
	# Interpolate the contour segments
	c_interps = []
	for x_seg, y_seg in zip(domains, ranges):
		c_interp = interp1d(x_seg, y_seg, 'linear')
		c_interps.append(c_interp)
	
	if peak_pos[0] > ir_pos[0]:
		tail = np.vstack((extrap_neg, ir_pos, peak_pos, extrap_pos))
	else:
		tail = np.vstack((extrap_pos, ir_pos, peak_pos, extrap_neg))
	
	tail_interp = interp1d(tail.T[0], tail.T[1], 'linear')
	
	# Find the intersections of the contours and tail
	x_intersects, y_intersects = [], []
	for ix, c_interp in enumerate(c_interps):
		x_intersect = curve_intersect(tail_interp, c_interp, domains[ix][0], domains[ix][-1])
		y_intersect = c_interp(x_intersect)
		x_intersects.append(x_intersect)
		y_intersects.append(y_intersect)
	
	intersects = np.vstack((np.hstack(x_intersects), np.hstack(y_intersects))).T
	
	# Return the maximum separation between host and edge
	intersects_sky = coord.SkyCoord(w.wcs_pix2world(intersects,1), unit=(u.deg,u.deg), frame='icrs')
	return max(ir.separation(intersects_sky))

def curve_intersect(fun1, fun2, xmin, xmax):
	'''
	Finds the intersection of two curves, bounded in [xmin, xmax]
	Returns an array of x values
	'''
	diff = lambda x: fun1(x)-fun2(x)
	x_range = np.linspace(xmin, xmax, 100)
	m_sign = np.sign(diff(x_range)).astype(int)
	roots = x_range[np.where(m_sign[1:] - m_sign[:-1] != 0)[0] + 1]
	
	# If they don't cross, return None
	if len(roots) == 0:
		return np.array([])
	
	# If they cross exactly once, find the global solution
	elif len(roots) == 1:
		return np.array([brentq(diff, xmin, xmax)])
	
	# If they cross multiple times, find the local solution between each root
	else:
		limits = np.concatenate(([xmin], roots, [xmax]))
		intersections = np.empty(len(limits)-2)
		for ix in range(len(intersections)):
			intersections[ix] = brentq(diff, limits[ix], limits[ix+1])
		return intersections

def peak_edge_ratio(w, ir, peaks, tails):
	'''
	Calculate the ratio of the distance to the peak and to the edge of each tail (measured on the sky)
	'''
	ratios = []
	for peak, tail in zip(peaks, tails):
		peak_pos = coord.SkyCoord(peak['ra'], peak['dec'], unit=(u.deg,u.deg), frame=('icrs'))
		ratios.append(ir.separation(peak_pos).deg/tail.deg)
	return ratios