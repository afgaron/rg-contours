# Radio galaxy contours

This repository contains a collection of algorithms for analyzing radio image contour data, based on code originally written for the Radio Galaxy Zoo pipeline. Two implementations of storing contour data as tree data structures are included:

- [`Contour_node`](contour_node.py) is a specialized tree that takes as input a contour file and a FITS image, and calculates physical properties for each island of emission (referred to as a "component") within the radio galaxy. These properties include integrated flux, integrated flux error, and the location of each emission peak (defined as the brightest pixel within each innermost contour). For integrated flux measurements, the beam size needs to be manually updated within the source code.

- [`Contour_path_object`](contour_path_object.py) is a more general tree that only requires a contour file and a host galaxy position. It does not store derived quantities for the radio galaxy, but contains functions used in calculating the bending angle, orientation angle, and tail lengths of the galaxy.

All of the functions used for processing the trees and deriving useful quantities are stored in [`processing`](processing.py), with a few miscellaneous helper functions used in both building the trees and processing the results stored in [`misc_functions`](misc_functions.py). Instructions on calling the two user-facing functions, `get_radio()` and `get_bending()`, are listed in [`example_usage`](example_usage.py).

## Inputs and outputs

An example FITS image of a radio galaxy and its corresponding contour file, stored as a JSON, are provided. These data correspond to the Radio Galaxy Zoo source [displayed here](https://radiotalk.galaxyzoo.org/#/subjects/ARG00026qx). The output of this code is a dictionary; an annotated print-out for the example data is provided below. Additional discussion about the definitions used for bending-related quantities can be found in [Garon et al. (2019)](https://iopscience.iop.org/article/10.3847/1538-3881/aaff62).

`{'bending': {'morphology': 'triple',`

The morphology of the source is based on the number of peaks in the radio emission: double sources have two peaks and triple sources have three. The bending algorithm can only take these two options as input right now.

`             'using_contour': {'asymmetry': 1.091938058857924,`

Two methods are implemented for determining the tail lengths and directions of a radio galaxy. The first is the "contour" method, in which a vector is drawn from the host galaxy position to the most distant point on the contours in either lobe. When using this method, outer contours are dropped from the data until there are two disjoint components remaining, neither of which contains the location of the host galaxy. This has the downside of sometimes adding a higher signal-to-noise threshold to the radio data than was used in generating the contours, and may even apply different thresholds to the two tails.

Asymmetry is defined as the ratio of the length of the longer tail to that of the shorter tail. It is always greater than or equal to 1.

`                               'bending_angle': 2.5097559521478154,`

The deviation of the two tails from co-linearity, measured in degrees (e.g., a perfectly straight source has a bending angle of 0°).

`                               'bisector': 101.53641087015717,`

The position angle of the angle bisector of the bending angle, measured in degrees east of north. This quantity can be used to describe the orientation on the sky of the radio galaxy.

`                               'pos_angle_0': 12.791288846231073,`
`                               'pos_angle_1': 190.28153289408326,`

Position angle of the two tails.

`                               'ratio_0': 0.7300105241522801,`
`                               'ratio_1': 0.7971262746886963,`

Ratio of the separation between the peak position and host position to the length of the tail. This corresponds to the Fanaroff-Riley classification, where a ratio of less than 0.5 is an FR-I-type tail and a ratio of greater than 0.5 is an FR-II-type tail.

`                               'size_deg': 0.015225278416527347,`

Sum of the lengths of the two tails, measured in degrees.

`                               'tail_deg_0': 0.006651388043754767,`
`                               'tail_deg_1': 0.00857389037277258},`

Length of the two tails, measured in degrees.

`             'using_peaks': {'asymmetry': 1.0928916153383081,`

The second method for determining the tail lengths and directions is the "peak" method, in which a vector is drawn from the host galaxy position to the radio peak in either lobe. Dropping outer contours like in the "contour" method has no effect on the morphology measured this way. All other quantities are calculated the same way as in the "contour" method, but using these new tail vectors.

`                             'bending_angle': 3.331708314934486,`
`                             'bisector': 101.38986224002792,`
`                             'pos_angle_0': 13.055716397495148,`
`                             'pos_angle_1': 189.72400808256066,`
`                             'ratio_0': 0.7300934243602069,`
`                             'ratio_1': 0.7979129818969033,`
`                             'size_deg': 0.015216069688823034,`
`                             'tail_deg_0': 0.006650632795955744,`
`                             'tail_deg_1': 0.00856543689286729}},`

` 'contour_file': '52af7d53eb9a9b05ef000809.json',`
` 'fits_file': 'FIRSTJ124646.7+224554.fits',`
` 'host_data': {'dec': 22.7615058, 'ra': 191.6946556},`

The contour file, FITS file, and host position are provided by the user.

` 'radio': {'components': [{'angular_extent': 32.36102391524036,`

Each radio component is located within a bounding box, a rectangle that circumscribes the radio contours. Angular extent is the length of the diagonal across that rectangle, measured in arcsec.

`                           'dec_range': [22.760603460952698,`
`                                         22.768125754011891],`

RA range and Dec range are the limits of the bounding box in RA and Dec.

`                           'flux': 36.556931587803035,`
`                           'flux_err': 0.28594502227989665,`

The integrated flux and corresponding uncertainty for the pixels contained within the contour, measured in mJy.

`                           'ra_range': [191.69274645933797,`
`                                        191.69808361969939],`
`                           'solid_angle': 268.46874254441764},

The projected area contained within the contours, measured in square arcsec.
`
`                          {'angular_extent': 27.22140162075618,`
`                           'dec_range': [22.752938330114109,`
`                                         22.759563725543924],`
`                           'flux': 34.352130434458338,`
`                           'flux_err': 0.26395570541452279,`
`                           'ra_range': [191.69142736631946,`
`                                        191.69537921470882],`
`                           'solid_angle': 228.76561864700378}],`

`           'dec': 22.760532042062998,`

RA and Dec for the radio galaxy refer to the center of the bounding box containing all emission.

`           'max_angular_extent': 32.104074813975352,`

The length of the diagonal across the bounding box containing all emission, measured in arcsec.

`           'number_components': 2,`

The number of isolated regions of radio emission in the galaxy.

`           'number_peaks': 3,`

The number of local maxima in the radio emission, maximum one per innermost contour level.

`           'outermost_level': 0.846972844901181,`

The flux limit used for the outermost contour level, measured in mJy/beam.

`           'peak_flux_err': 0.02399597321950207,`

The uncertainty in peak flux measurements, measured in mJy/beam

`           'peaks': [{'dec': 22.76623586618534,`
`                      'flux': 16.62060059607029,`
`                      'ra': 191.6958451445872},`

For each radio peak, its location in RA and Dec and its value in mJy/beam.

`                     {'dec': 22.761648521120179,`
`                      'flux': 1.6361365560442209,`
`                      'ra': 191.69462003720014},`
`                     {'dec': 22.754769514196084,`
`                      'flux': 23.802017793059349,`
`                      'ra': 191.69340381298568}],`

`           'ra': 191.69475549300944,`
`           'total_flux': 70.909062022261367,`
`           'total_flux_err': 0.38914929035979101,`

The total flux is the sum of the fluxes of the components, and the uncertainty is the sum in quadrature of the component flux uncertainties.

`           'total_solid_angle': 497.23436119142139}}`

The total solid angle is the sum of the component solid angles.

## Additional help

Background on the algorithms used here can be found in [Banfield et al. (2015)](https://academic.oup.com/mnras/article/453/3/2326/1075547) and [Garon et al. (2019)](https://iopscience.iop.org/article/10.3847/1538-3881/aaff62). For additional help, please contact Avery Garon at [garo0040@umn.edu](mailto:garo0040@umn).
