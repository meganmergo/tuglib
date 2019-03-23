#!/usr/bin/env python

__all__ = ['FitsCollection', 'image_combine', 'bias_combine',
           'dark_combine', 'flat_combine']

from .helper import FitsCollection
from ccdproc import CCDData


# Generic image combine function.
# This will be bases for all bias/dark/flat combine methods
def image_combine(images, method='median', masks=None,
                  trim=None, return_ccddata=True,
                  output_fname='master_image.fits',
                  gain=None, read_noise=None, **kwargs):

    pass


def bias_combine():
    pass


def dark_combine():
    pass


def flat_combine():
    pass


