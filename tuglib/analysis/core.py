#!/usr/bin/env python

__all__ = ['phot']

import types

from photutils import aperture_photometry
from photutils import CircularAperture, CircularAnnulus
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground
from photutils.utils import calc_total_error
from astropy.table import vstack
from astropy.nddata import StdDevUncertainty
from astropy import units as u

import numpy as np


# Generic photometry function
def phot(images, positions=None, aperture_radius=6.0, annulus=8.0, dannulus=10.0, method='exact', bg_method='median',
         bg_box_size=50,
         gain=1.1,
         output=None):
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


    aperture_radius : list or float
        The aperture radius of a source.

        Default is '3.0'.

    annulus : float
        The circular inner annulus aperture radius of a source.

        Default is '6.0'.

    dannulus : float
        The circular outer annulus aperture radius of a source.

        Default is '8.0'.

    method : {'exact', 'center', 'subpixel'}, optional
        The method used to determine the overlap of the aperture on the
        pixel grid.  Not all options are available for all aperture
        types.  Note that the more precise methods are generally slower.
        The following methods are available:

            * ``'exact'`` (default):
                The the exact fractional overlap of the aperture and
                each pixel is calculated.  The returned mask will
                contain values between 0 and 1.

            * ``'center'``:
                A pixel is considered to be entirely in or out of the
                aperture depending on whether its center is in or out of
                the aperture.  The returned mask will contain values
                only of 0 (out) and 1 (in).

            * ``'subpixel'``:
                A pixel is divided into subpixels (see the ``subpixels``
                keyword), each of which are considered to be entirely in
                or out of the aperture depending on whether its center
                is in or out of the aperture.  If ``subpixels=1``, this
                method is equivalent to ``'center'``.  The returned mask
                will contain values between 0 and 1.

    bg_method: {'mean', 'median'}, optional
        The statistic used to calculate the background.
        All measurements are sigma clipped.

    bg_box_size: int
    Selecting the box size requires some care by the user. The box size
    should generally be larger than the typical size of sources in the
    image, but small enough to encapsulate any background variations. For
    best results, the box size should also be chosen so that the data are
    covered by an integer number of boxes in both dimensions. More information:
    https://github.com/Gabriel-p/photpy/blob/master/IRAF_compare/aperphot.py

    gain: float
        CCD gain value.

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
    >>> phot(fits_images, positions=[(1400.65, 1545.78)], aperture_radius=10.0, annulus=12.0, dannulus=15)
    """

    if not isinstance(images, (list, types.GeneratorType)):
        raise TypeError(
            "'images' should be a 'ccdproc.CCDData' object.")

    if positions is not None:
        if not isinstance(positions, (list, type(None))):
            raise TypeError(
                "'positions' should be 'list' or 'None' object.")

    if not isinstance(aperture_radius, (list, float, int)):
        raise TypeError("'radius' should be a 'list' or 'float' object.")

    if not isinstance(annulus, (float, int)):
        raise TypeError("'annulus' should be a 'int' or 'float' object.")

    if not isinstance(dannulus, (float, int)):
        raise TypeError("'dannulus' should be a 'int' or 'float' object.")

    if not isinstance(method, str):
        raise TypeError("'method' should be a 'str' object.")

    if not isinstance(bg_method, str):
        raise TypeError("'bg_method' should be a 'str' object.")

    if not isinstance(bg_box_size, int):
        raise TypeError("'bg_box_size' should be a 'int' or 'float' object.")

    if not isinstance(gain, (float, int)):
        raise TypeError("'gain' should be a 'float' object.")

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

        box_xy = (bg_box_size, bg_box_size)
        sigma_clip = SigmaClip(sigma=3.)
        bkg_estimator = MedianBackground()
        bkg = Background2D(ccd, box_xy, filter_size=(3, 3),
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

        ccd.uncertainty = StdDevUncertainty(calc_total_error(ccd.data, bkg.background, gain))

        aperture = CircularAperture(positions, r=aperture_radius)
        annulus_aperture = CircularAnnulus(positions, r_in=annulus, r_out=dannulus)
        annulus_masks = annulus_aperture.to_mask(method=method)

        bkg_phot = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(ccd)
            annulus_data_1d = annulus_data[mask.data > 0]
            mean_sigclip, median_sigclip, stddev_sigclip = sigma_clipped_stats(annulus_data_1d)

            if bg_method == "mean":
                sigclip = mean_sigclip
            elif bg_method == "median":
                sigclip = median_sigclip
            else:
                sigclip = median_sigclip

            bkg_phot.append(sigclip)

        bkg_phot = np.array(sigclip)
        phot_table = aperture_photometry(ccd, aperture)
        phot_table['annulus_{}'.format(bg_method)] = bkg_phot * u.adu
        phot_table['aper_bkg'] = bkg_phot * aperture.area * u.adu
        phot_table['residual_aperture_sum'] = phot_table['aperture_sum'] - phot_table['aper_bkg']

        phot_table['flux'] = phot_table['residual_aperture_sum'] / float(exptime)
        phot_table['flux'].info.format = '%.3f'
        phot_table['JD'] = float(ccd.meta['JD'])
        phot_table['JD'].info.format = '%.8f'
        phot_table['JD'].unit = 'd'
        phot_table['aperture_sum'].info.format = '%.3f'
        phot_table['aperture_sum_err'].info.format = '%.3f'
        phot_table['annulus_{}'.format(bg_method)].info.format = '%.3f'
        phot_table['aper_bkg'].info.format = '%.3f'
        phot_table['residual_aperture_sum'].info.format = '%.3f'
        phot_table = flux2mag(phot_table)

        phot_table['mag_err'] = 1.0857 * (phot_table['aperture_sum_err'] / phot_table['residual_aperture_sum']) * u.mag
        phot_table['mag_err'].info.format = '%.3f'
        ccds_phot_list.append(phot_table)

    return vstack(ccds_phot_list)


def flux2mag(phot_table, zmag=25.):
    phot_table['mag'] = (zmag - 2.5 * np.log10(np.array(phot_table['flux']))) * u.mag
    phot_table['mag'].info.format = '%.3f'
    return phot_table
