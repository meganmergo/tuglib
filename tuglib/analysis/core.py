#!/usr/bin/env python

__all__ = ['phot', 'test_phot']

import types

import numpy as np

from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry

from astropy.table import Table, vstack

from ccdproc import CCDData

from .helper import radec_to_pixel_coordinate, get_object_coordinate_from_simbad


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


# Use at your own risk! Experimental :)
def test_phot(imgs, sources, parameters):
    """Perform aperture photometry.

    Parameters
    ----------
    imgs : list
        List of 'fits' image file names.

    sources : list
        List of astronomical object names.

    parameters : tuple
        Aperture photmetry parameters (r, r_in, r_out).
        r: inner aperture
        r_in: sky inward radius
        r_out: sky outward radius

    Examples
    --------
    >>> imgs = sorted(glob('*.fits'))
    >>> r, r_in, r_out = 9, 12, 15
    >>> parameters = (r, r_in, r_out)
    >>> sources = ['V2477 Cyg', 'TYC 3945-1069-1', 'TYC 3945-1153-1',\
                   '1SWASP J201900.18+563607.3', 'TYC 3945-1162-1',\
                   'TYC 3945-1197-1']
    >>>
    >>> r = test_phot(imgs, sources, parameters)
    >>> diff_mags = -2.5 * np.log10(r['V2477 Cyg'] - r['TYC 3945-1197-1']) + 30
    >>>
    >>> plt.grid()
    >>> plt.title('Light Curve')
    >>> plt.xlabel('JD')
    >>> plt.ylabel('Differential Magnitude')
    >>> plt.gca().invert_yaxis()
    >>> plt.plot(r['JD'], diff_mags, 'ro')
    >>> plt.show()
    """

    fluxes = list()
    sky_coords = list()

    for source in sources:
        sky_coords.append(get_object_coordinate_from_simbad(source))

    r, r_in, r_out = parameters

    for img in imgs:
        ccd = CCDData.read(img)
        exptime = ccd.meta['EXPTIME']

        positions = list()
        for coord in sky_coords:
            x, y = radec_to_pixel_coordinate(ccd, coord)
            positions.append((x, y))

        circular_apertures = CircularAperture(positions, r=r)
        annulus_apertures = CircularAnnulus(positions, r_in=r_in, r_out=r_out)

        apers = [circular_apertures, annulus_apertures]
        phot_table = aperture_photometry(ccd.data, apers)

        bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area()
        bkg_sum = bkg_mean * circular_apertures.area()

        final_sum = phot_table['aperture_sum_0'] - bkg_sum
        phot_table['residual_aperture_sum'] = final_sum
        phot_table['residual_aperture_sum'].info.format = '%.8g'

        flux = phot_table['residual_aperture_sum']
        flux = flux.data / float(exptime)

        flux = np.insert(flux, 0, float(ccd.meta['JD']))
        fluxes.append(flux)

    names = ['JD'] + sources
    result = Table(rows=fluxes, names=names)

    return result
