#!/usr/bin/env python

__all__ = ['phot']

import types

from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import CircularAperture, CircularAnnulus
from astropy.table import Table, vstack
from astropy import units as u

import numpy as np

from .helper import radec_to_pixel_coordinate, get_object_coordinate_from_simbad


# Generic photometry function
def phot(images, positions=None, radius=6.0, annulus=8.0, dannulus=10.0, output=None):
    """
    Basic circular aperture photometry of given images.

    Parameters
    ----------
    images : generator or list of 'ccdproc.CCDData'
        Images to be combined.

    positions : list
        The positions should be either a single tuple of (x, y), a list of (x, y) tuples, or
        an array with shape Nx2, where N is the number of positions.

        Default is 'None'.


    radius : list or float
        The aperture radius of a source.

        Default is '3.0'.

    annulus : float
        The circular inner annulus aperture radius of a source.

        Default is '6.0'.

    dannulus : float
        The circular outer annulus aperture radius of a source.

        Default is '8.0'.

    output : None or str
        If it is None, function returns just Table.
        If it is 'str', function creates photometry result file.


    Returns
    -------
    'A table object or file'

    Examples
    --------
    >>> from tuglib.io import FitsCollection
    >>> from tuglib.analysis import phot
    >>> location = '/Users/ykilic/Desktop/data/'
    >>> images = FitsCollection(location, file_extension="fit")
    >>> fits_images = images.ccds(IMAGETYP='object')
    >>> phot(fits_images, positions=[(1400.65, 1545.78)], radius=6.0, annulus=8.0, dannulus=10)
    """

    if not isinstance(images, (list, types.GeneratorType)):
        raise TypeError(
            "'images' should be a 'ccdproc.CCDData' object.")

    if positions is not None:
        if not isinstance(positions, (list, type(None))):
            raise TypeError(
                "'positions' should be 'list' or 'None' object.")

    if not isinstance(radius, (list, float, int)):
        raise TypeError("'radius' should be a 'list' or 'float' object.")

    if not isinstance(annulus, (float, int)):
        raise TypeError("'annulus' should be a 'int' or 'float' object.")

    if not isinstance(dannulus, (float, int)):
        raise TypeError("'dannulus' should be a 'int' or 'float' object.")

    if not isinstance(output, (type(None), type(str))):
        raise TypeError("'output' should be 'None' or 'str' objects.")

    if isinstance(images, types.GeneratorType):
        try:
            ccds = list(images)
        except IndexError:
            return None
    else:
        ccds = images

    ccds_phot_list = list()

    for ccd in ccds:
        exptime = ccd.meta['EXPTIME']

        aperture = CircularAperture(positions, r=radius)
        annulus_aperture = CircularAnnulus(positions, r_in=annulus, r_out=dannulus)
        annulus_masks = annulus_aperture.to_mask(method='center')

        bkg_median = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(ccd)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        bkg_median = np.array(bkg_median)
        phot_table = aperture_photometry(ccd, aperture)
        phot_table['annulus_median'] = bkg_median * u.adu
        phot_table['aper_bkg'] = bkg_median * aperture.area * u.adu
        phot_table['aper_sum_bkgsub'] = phot_table['aperture_sum'] - phot_table['aper_bkg']

        phot_table['flux'] = phot_table['aper_sum_bkgsub'] / float(exptime)
        phot_table['JD'] = float(ccd.meta['JD'])
        phot_table['JD'].info.format = '%.8f'
        phot_table['JD'].unit = 'd'
        phot_table['aperture_sum'].info.format = '%.3f'
        phot_table['annulus_median'].info.format = '%.3f'
        phot_table['aper_bkg'].info.format = '%.3f'
        phot_table['aper_sum_bkgsub'].info.format = '%.3f'
        phot_table['flux'].info.format = '%.3f'

        ccds_phot_list.append(phot_table)

    return vstack(ccds_phot_list)