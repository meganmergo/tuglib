#!/usr/bin/env python

__all__ = ['FitsCollection']

from os import path
from glob import glob

import numpy as np

from astropy import units as u
from astropy.table import Table
from astropy.io import fits

from ccdproc import ImageFileCollection, CCDData


FILE_EXTENSIONS = ('fit', 'fits', 'fit.gz', 'fits.gz', 'fit.zip', 'fits.zip')


# Another helper method for FitsCollection class.
def get_fits_header(filename):
    """Transposes fits header column by column.

    Parameters
    ----------
    filename : str
        Fits file to be trasposed.

    Returns
    -------
    result : list of list
        [keywords, values, comments, dtypes]
    """

    h = fits.getheader(filename)

    keywords = list(h.keys())
    values = list(h.values())
    comments = list(h.comments)

    dtypes = list()
    for x in values:
        if isinstance(x, bool):
            dtypes.append('bool')
        elif isinstance(x, int):
            dtypes.append('int')
        elif isinstance(x, float):
            dtypes.append('float')
        else:
            dtypes.append('U32')

    return keywords, values, comments, dtypes


# It doesn't work well!
# I will write completely new class for this purpose.
# Helper class for ccdproc.ImageFileCollection
class FitsCollectionTest(ImageFileCollection):

    def __init__(self, **kwargs):
        if 'location' in kwargs:
            directory = path.join(kwargs['location'], '**', '*.*')

            filenames = glob(directory, recursive=True)
            kwargs['filenames'] = filenames

        super(FitsCollectionTest, self).__init__(**kwargs)


# New class for recursive image collection.
class FitsCollection(object):
    """FITS Image Collection
    (A ccdproc.ImageFileCollection alternative)

    It performs recursive fits image search in given directory.
    """

    def __init__(self, location, file_extension=None,
                 start_init=True, **kwargs):

        if not isinstance(location, str):
            raise TypeError(
                "'location' should be a 'str' object.")

        if not isinstance(file_extension, (type(None), str, tuple, list)):
            raise TypeError(
                "'file_extension' should be "
                "'None', 'str', 'tuple' or 'list' object.")

        if file_extension not in FILE_EXTENSIONS:
            raise ValueError(
                "'file_extension' should be "
                "'fit, fits, fit.gz, fits.gz, fit.zip, fits.zip' type.")

        if not isinstance(start_init, bool):
            raise TypeError("'start_init' should be a 'bool' object.")

        self.location = location
        self._file_extension = file_extension
        self.filenames = list()
        self._start_init = start_init

        self.collection = None

        if location is not None:
            if self._start_init:
                self.init()

    def init(self):
        directory = list()

        if self._file_extension is None:
            for ext in FILE_EXTENSIONS:
                directory += path.join(self.location, '**', '*.' + ext)
        elif isinstance(self._file_extension, (tuple, list)):
            for ext in self._file_extension:
                directory += path.join(self.location, '**', '*.' + ext)
        else:  # str
            directory = path.join(
                self.location, '**', '*.' + self._file_extension)

        self.filenames = sorted(glob(directory, recursive=True))

        r = get_fits_header(self.filenames[0])
        names, dtypes = r[0], r[3]
        names = ['filename'] + names
        dtypes = ['U64'] + dtypes

        self.collection = Table(names=names, dtype=dtypes)

        for filename in self.filenames:
            r = get_fits_header(filename)
            values = [filename] + r[1]

            self.collection.add_row(values)  # values

    def __str__(self):
        pass

    def __repr__(self):
        pass
