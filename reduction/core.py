#!/usr/bin/env python

__all__ = ['image_combine', 'bias_combine', 'dark_combine', 'flat_combine']

from astropy.stats import mad_std
from ccdproc import ImageFileCollection
from ccdproc import combine
import numpy as np

# Generic image combine function.
# This will be bases for all bias/dark/flat combine methods


def image_combine(img_list, method='average', masks=None,
                  trim=None, return_ccddata=True,
                  output_fname='master_image.fits',
                  gain=None, read_noise=None, **kwargs):

    combined_images = combine(img_list,
                              output_file=output_fname,
                              method=method,
                              sigma_clip=True,
                              sigma_clip_low_thresh=5,
                              sigma_clip_high_thresh=5,
                              sigma_clip_func=np.ma.median,
                              signma_clip_dev_func=mad_std)

    if isinstance(output_fname, str) or isinstance(output_fname, unicode):
        return True
    elif output_fname is None:
        return combined_images


def bias_combine():
    pass


def dark_combine():
    pass


def flat_combine():
    pass
