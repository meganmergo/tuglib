#!/usr/bin/env python

__all__ = ['FitsCollection', 'make_mask', 'image_combine', 'bias_combine',
           'dark_combine', 'flat_combine']

import numpy as np

from .helper import FitsCollection, make_mask
from ccdproc import CCDData, trim_image, combine


# Generic image combine function.
# This will be bases for all bias/dark/flat combine methods
def image_combine(images, method='median', output='master_image.fits',
                  masks=None, trim=None, gain=None, read_noise=None,
                  return_ccddata=True, **kwargs):

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

    ccds = images.ccds(**kwargs)

    if masks is not None:
        shape = ccds[0].data.shape
        mask = make_mask(shape, masks[0])
        for m in masks[1:]:
            tmp_mask = make_mask(shape, m)
            mask |= tmp_mask
    else:
        mask = None

    tmp_ccds = list()
    if trim is not None:
        for ccd in ccds:
            ccd = trim_image(ccd, trim)
            tmp_ccds.append(ccd)





def bias_combine():
    pass


def dark_combine():
    pass


def flat_combine():
    pass


