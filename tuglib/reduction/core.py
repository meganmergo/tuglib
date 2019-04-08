#!/usr/bin/env python

__all__ = ['FitsCollection', 'make_mask', 'image_combine', 'bias_combine',
           'dark_combine', 'flat_combine']

import types

from .helper import FitsCollection, make_mask
from ccdproc import CCDData, trim_image, combine, subtract_bias


# Generic image combine function.
# This will be bases for all bias/dark/flat combine methods
def image_combine(images, method='median', output=None, masks=None, trim=None):
    """
    Generic Image Combine Method.

    Parameters
    ----------
    images : generator or list of 'ccdproc.CCDData'
        Images to be combined.

    method : str
        Method to combine images:

        - average : To combine by calculating the average.
        - median : To combine by calculating the median.
        - sum : To combine by calculating the sum.

        Default is 'median'.

    output : None or str
        If it is None, function returns just 'ccdproc.CCDData'.
        If it is 'str', function returns 'ccdproc.CCDData' and creates file.

    masks : str, list of str or optional
        Area to be masked.

    trim : str or optional
        Trim section.

    Returns
    -------
    'ccdproc.CCDData'
        Combined Images.

    Examples
    --------

    >>> from tuglib.reduction import FitsCollection, image_combine
    >>>
    >>> path = '/home/user/data/'
    >>> masks = ['[:, 1023:1025]', '[:1023, 56:58]']
    >>> trim = '[:, 24:2023]'
    >>>
    >>> images = FitsCollection(location=path, gain=0.57, read_noise=4.11)
    >>> ccds = images.ccds(OBJECT='Star', masks=masks, trim=trim)
    >>>
    >>> master_image = combine(ccds, method='median', output='master.fits')
    """

    if not isinstance(images, (list, types.GeneratorType)):
        raise TypeError(
            "'images' should be a 'ccdproc.CCDData' object.")

    if method not in ('average', 'median', 'sum'):
        raise ValueError(
            "'method' should be 'average', 'median' or 'sum'.")

    if masks is not None:
        if not isinstance(masks, (str, list, type(None))):
            raise TypeError(
                "'masks' should be 'str', 'list' or 'None' object.")

    if trim is not None:
        if not isinstance(trim, str):
            raise TypeError("'trim' should be a 'str' object.")

    if not isinstance(output, (type(None), str)):
        raise TypeError("'output' should be 'None' or 'str' objects.")

    if isinstance(images, types.GeneratorType):
        ccds = list(images)
    else:
        ccds = images

    mask = None
    if masks is not None:
        shape = ccds[0].shape
        mask = make_mask(shape, masks)

    for i, ccd in enumerate(ccds):
        if mask is not None:
            ccd.mask = mask

        ccd = trim_image(ccd, trim)

        ccds[i] = ccd

    master_ccd = combine(ccds, method=method)

    if output is not None:
        master_ccd.write(output, overwrite=True, output_verify='ignore')

    return master_ccd


def bias_combine(images, method='median', output=None,
                 masks=None, trim=None):
    """
    Bias Combine.

    Parameters
    ----------
    images : generator or list of 'ccdproc.CCDData'
        Images to be combined.

    method : str
        Method to combine images:

        - average : To combine by calculating the average.
        - median : To combine by calculating the median.
        - sum : To combine by calculating the sum.

        Default is 'median'.

    output : None or str
        If it is None, function returns just 'ccdproc.CCDData'.
        If it is 'str', function returns 'ccdproc.CCDData' and creates file.

    masks : str, list of str or optional
        Area to be masked.

    trim : str or optional
        Trim section.

    Returns
    -------
    master_ccd : ccdproc.CCDData
        Combined Images.

    Examples
    --------

    >>> from tuglib.reduction import FitsCollection, bias_combine
    >>>
    >>> path = '/home/user/data/'
    >>> masks = ['[:, 1023:1025]', '[:1023, 56:58]']
    >>> trim = '[:, 24:2023]'
    >>>
    >>> images = FitsCollection(location=path, gain=0.57, read_noise=4.11))
    >>> bias_ccds = images.ccds(OBJECT='BIAS', trim=trim, masks=masks)
    >>>
    >>> master_bias = bias_combine(bias_ccds, method='median',
                                   output='master_bias.fits')
    """

    return image_combine(images, method=method, output=output,
                         masks=masks, trim=trim)


def dark_combine(images, master_bias=None, method='median',
                 output=None, masks=None, trim=None):

    if not isinstance(images, (list, types.GeneratorType)):
        raise TypeError(
            "'images' should be a 'ccdproc.CCDData' object.")

    if not isinstance(master_bias, CCDData):
        raise TypeError(
            "'master_bias' should be a 'ccdproc.CCDData' object.")

    if method not in ('average', 'median', 'sum'):
        raise ValueError(
            "'method' should be 'average', 'median' or 'sum'.")

    if masks is not None:
        if not isinstance(masks, (str, list, type(None))):
            raise TypeError(
                "'masks' should be 'str', 'list' or 'None' object.")

    if trim is not None:
        if not isinstance(trim, str):
            raise TypeError("'trim' should be a 'str' object.")

    if not isinstance(output, (type(None), str)):
        raise TypeError("'output' should be 'None' or 'str' objects.")

    if isinstance(images, types.GeneratorType):
        ccds = list(images)
    else:
        ccds = images

    mask = None
    if masks is not None:
        shape = ccds[0].shape
        mask = make_mask(shape, masks)

    for i, ccd in enumerate(ccds):
        if mask is not None:
            ccd.mask = mask

        ccd = trim_image(ccd, trim)

        ccds[i] = ccd

    raw_darks = dict()
    for ccd in ccds:
        exp_time = ccd.meta['EXPTIME']

        if raw_darks.get(exp_time) is None:
            raw_darks[exp_time] = list()
            raw_darks[exp_time].append(ccd)
        else:
            raw_darks[exp_time].append(ccd)

    master_darks = list()
    for darks in raw_darks.values():
        bias_subtracked_darks = list()
        for dark in darks:
            bias_subtracked_darks.append(subtract_bias(dark, master_bias))

        master_darks.append(image_combine(bias_subtracked_darks, method=method))

    if output is not None:
        for master_dark in master_darks:
            filename = output + '_' + str(master_dark.meta['EXPTIME']) + '.fits'
            master_dark.write(filename, overwrite=True, output_verify='ignore')

    return master_darks


def flat_combine():
    pass
