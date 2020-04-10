#!/usr/bin/env python

__all__ = ['radec_to_pixel_coordinate', 'get_object_coordinate_from_simbad']

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord

from astroquery.simbad import Simbad


def radec_to_pixel_coordinate(ccd, coord):
    """Convert RA and Dec to x and y pixel coordinate.
    'ccd' must have 'wcs' keywords.

    Parameters
    ----------

    ccd : CCDData
        CCDData object.

    coord : tuple
        Contain RA and Dec coordinates.

    Returns
    -------

    pos : np.array
        x and y coordinates.

    Examples
    --------

    >>> coord = ('13:12:11.10', '+36:49:27.13')
    >>> p = radec_to_pixel_coordinate(ccd, coord)
    """

    p = np.array([[coord.ra.value[0], coord.dec.value[0]]])
    pos = ccd.wcs.wcs_world2pix(p, 1)[0]
    return pos


def get_object_coordinate_from_simbad(name):
    """Get coordinate of object using simbad service.

    Parameters
    ----------

    name : str
        Name of object.

    Returns
    -------

    coord : tuple
        RA and Dec of object.

    Examples
    --------

    >>> coord = get_object_coordinate_from_simbad('W UMa')
    """
    obj = Simbad.query_object(name)
    coord = SkyCoord(obj['RA'], obj['DEC'], unit=(u.hourangle, u.deg))
    return coord
