#!/usr/bin/env python

__all__ = ['FitsCollection', 'make_mask', 'image_combine', 'bias_combine',
           'dark_combine', 'flat_combine']

from .helper import FitsCollection, make_mask
from ccdproc import combine


# Generic image combine function.
# This will be bases for all bias/dark/flat combine methods
def image_combine(images, method='median', output=None,
                  masks=None, trim=None, **kwargs):
    """
    Generic Image Combine Method.

    Parameters
    ----------
    images : FitsCollection
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

    **kwargs: dict
        Additional keywords are used to filter FitCollection.

    Returns
    -------
    master_ccd : ccdproc.CCDData
        Combined Images.

    Examples
    --------

    >>> from tuglib.reduction import FitsCollection, image_combine
    >>>
    >>> path = '/home/user/data/'
    >>> masks = ['[:, 1023:1025]', '[:1023, 56:58]']
    >>> trim = '[:, 24:2023]'
    >>>
    >>> images = FitsCollection(location=path, gain=0.57, read_noise=4.11))
    >>> master_image = image_combine(images, method='median',
                                     output='master.fits',
                                     masks=masks, trim=trim,
                                     OBJECT='Star')
    """

    if not isinstance(images, FitsCollection):
        raise TypeError(
            "'images' should be a 'FitsCollection' object.")

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

    ccds = images.ccds(trim=trim, masks=masks, **kwargs)

    master_ccd = combine(ccds, method=method)

    if output is not None:
        master_ccd.write(output, overwrite=True, output_verify='ignore')

    return master_ccd


def bias_combine(bias_images, method='median', output=None,
                 masks=None, trim=None, **kwargs):
    """
    Bias Combine.

    Parameters
    ----------
    bias_images : FitsCollection
        Bias images to be combined.

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

    **kwargs: dict
        Additional keywords are used to filter FitCollection.

    Returns
    -------
    master_ccd : ccdproc.CCDData
        Combined Images.

    Examples
    --------

    >>> from tuglib.reduction import FitsCollection, image_combine
    >>>
    >>> path = '/home/user/data/'
    >>> masks = ['[:, 1023:1025]', '[:1023, 56:58]']
    >>> trim = '[:, 24:2023]'
    >>>
    >>> images = FitsCollection(location=path, gain=0.57, read_noise=4.11))
    >>> master_bias = image_combine(images, method='median',
                                    output='master_bias.fits',
                                    masks=masks, trim=trim,
                                    OBJECT='BIAS')
    """

    return image_combine(bias_images, method=method, output=output,
                         masks=masks, trim=trim, **kwargs)


def dark_combine():
    pass


def flat_combine():
    pass
