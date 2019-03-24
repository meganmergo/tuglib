#!/usr/bin/env python

__all__ = ['FitsCollection', 'make_mask', 'image_combine', 'bias_combine',
           'dark_combine', 'flat_combine']

from .helper import FitsCollection
from ccdproc import CCDData


def make_mask(shape, area):
    """Make a mask specific to the given area.

    Parameters
    ----------

    shape : tuple
        Size of mask (m, n)

    area : str
        The area to be masked.

    Returns
    -------
    m : np.array
        Masked area.

    Examples
    --------

    >>> m = make_mask((1024, 1024), '[100:200, 23:55]')
    """

    m = np.full(shape, False, dtype=bool)
    exec('m' + area + ' = True')

    return m


# Generic image combine function.
# This will be bases for all bias/dark/flat combine methods
def image_combine(images, method='median', masks=None,
                  trim=None, return_ccddata=True,
                  output='master_image.fits',
                  gain=None, read_noise=None, **kwargs):

    if not isinstance(images, FitsCollection):
        raise TypeError(
            "'images' should be a 'FitsCollection' object.")

    if method not in ('average', 'median', 'sum'):
        raise ValueError(
            "'method' should be 'average', 'median' or 'sum'.")

    if masks is not None:
        if not isinstance(masks, list):
            raise TypeError("'masks' should be a 'list' object.")

    if trim is not None:
        if not isinstance(trim, str):
            raise TypeError("'trim' should be a 'str' object.")

    if not isinstance(return_ccddata, bool):
        raise TypeError("'return_ccddata' should be a 'bool' object.")

    if not isinstance(output, str):
        raise TypeError("'output' should be a 'str' object.")

    if not isinstance(gain, float):
        raise TypeError(
            "'gain' should be a 'float' object.")

    if not isinstance(read_noise, float):
        raise TypeError(
            "'read_noise' should be a 'float' object.")

    if masks is not None:
        shape = next(images.hdus(**kwargs)).shape
        mask = make_mask(shape, masks[0])
        for m in masks[1:]:
            tmp_mask = make_mask(shape, m)
            mask |= tmp_mask
    else:
        mask = None

    master_list = list()



def bias_combine():
    pass


def dark_combine():
    pass


def flat_combine():
    pass


