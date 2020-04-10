#!/usr/bin/env python

__all__ = ['phot']

import types

from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import CircularAperture, CircularAnnulus
from astropy.table import Table, vstack
from astropy import units as u

import numpy as np

# Generic photometry function
def phot(images, positions=None, radius=3.0, annulus=6.0, dannulus=8.0, output=None):
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
    >>>
    >>> path = '/home/user/data/'
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
        raise TypeError("'radius' should be a 'int' or 'float' object.")

    if not isinstance(dannulus, (float, int)):
        raise TypeError("'radius' should be a 'int' or 'float' object.")

    if not isinstance(output, (type(None), type(str))):
        raise TypeError("'output' should be 'None' or 'str' objects.")

    if isinstance(images, types.GeneratorType):
        try:
            ccds = list(images)
        except IndexError:
            return None
    else:
        ccds = images

    ccds_phot_list = []
    for ccd in ccds:
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
        aphot = aperture_photometry(ccd, aperture)
        aphot['annulus_median'] = bkg_median * u.adu
        aphot['aper_bkg'] = bkg_median * aperture.area * u.adu
        aphot['aper_sum_bkgsub'] = aphot['aperture_sum'] - aphot['aper_bkg']
        for col in aphot.colnames:
            aphot[col].info.format = '%.8g'  # for consistent table output
        ccds_phot_list.append(aphot)

    return vstack(ccds_phot_list)