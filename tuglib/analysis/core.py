#!/usr/bin/env python

__all__ = ['phot']

import types

from photutils import CircularAperture
from photutils import aperture_photometry
from astropy.table import Table, vstack

# Generic photometry function
def phot(images, positions=None, radius=3.0, output=None):
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
        phot_table = aperture_photometry(ccd, aperture)
        phot_table['aperture_sum'].info.format = '%.8g'
        ccds_phot_list.append(phot_table)

    return vstack(ccds_phot_list)