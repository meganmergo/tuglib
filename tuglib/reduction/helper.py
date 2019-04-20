#!/usr/bin/env python

__all__ = ['make_mask']

import numpy as np


def helper_mask(shape, area, fill=True):
    """
    Make a mask specific to the given area.

    Parameters
    ----------

    shape : tuple
        Size of mask (m, n)

    area : str
        The area to be masked.

    fill : bool, optional
        Masked area filled by 'fill' value. Default value 'True'.

    Returns
    -------
    m : np.array
        Masked area.

    Examples
    --------

    >>> m = helper_mask((1024, 1024), '[100:200, 23:55]')
    """

    if not isinstance(shape, tuple):
        raise TypeError(
            "'shape' shuld be a 'tuple' object.")

    if not isinstance(area, str):
        raise TypeError(
            "'masks' should be a 'str' object.")

    if not isinstance(fill, bool):
        raise TypeError("'fill' should be a 'bool' object.")

    m = np.full(shape, not fill, dtype=bool)

    if fill:
        exec('m' + area + ' = True')
    else:
        exec('m' + area + ' = False')

    return m


def make_mask(shape, masks, fill=True):
    """
    Make a mask specific to the given area.

    Parameters
    ----------

    shape : tuple
        Size of mask (m, n)

    masks : str or list of str
        The area to be masked.

    fill : bool, optional
        Masked area filled by 'fill' value. Default value 'True'.

    Returns
    -------
    m : np.array
        Masked area.

    Examples
    --------

    >>> shape = (1024, 1024)
    >>> masks = ['[100:200, 23:55]', '[511, :]']
    >>> m = make_mask(shape, masks)
    """

    if not isinstance(shape, tuple):
        raise TypeError(
            "'shape' shuld be a 'tuple' object.")

    if not isinstance(masks, (str, list)):
        raise TypeError(
            "'masks' should be 'str' or 'list' object.")

    if not isinstance(fill, bool):
        raise TypeError("'fill' should be a 'bool' object.")

    if isinstance(masks, list):
        mask = helper_mask(shape, masks[0], fill=fill)
        for m in masks[1:]:
            tmp_mask = helper_mask(shape, m, fill=fill)
            mask |= tmp_mask
    else:
        mask = helper_mask(shape, masks, fill=fill)

    return mask
