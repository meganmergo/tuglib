#!/usr/bin/env python

__all__ = ['FitsCollection']

import numpy as np

from os import path
from glob import glob

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
            dtypes.append('U64')

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

    Parameters
    ----------
    location : str
        Full path directory that include images.
        Example: '/home/user/data/20190325'

    file_extension : None, str, tuple or list
        Image extensions. Should be one of 'fit', 'fits', 'fit.gz',
        'fits.gz', 'fit.zip', 'fits.zip'

    Returns
    -------

    """

    def __init__(self, location, file_extension='fits'):

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

        self.location = location
        self._file_extension = file_extension
        self.filenames = list()

        self.collection = None

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

        db = dict()
        for filename in self.filenames:
            h = get_fits_header(filename)
            db[filename] = h

        cs = dict()
        for filename in self.filenames:
            for i, key in enumerate(db[filename][0]):
                cs[db[filename][0][i]] = db[filename][3][i]

        names = list(cs.keys())
        names = ['filename'] + names
        dtypes = list(cs.values())
        dtypes = ['U64'] + dtypes

        self.collection = Table(names=names, dtype=dtypes)

        for filename in self.filenames:
            self.collection.add_row()
            row_index = len(self.collection) - 1
            self.collection['filename'][row_index] = filename

            h = db[filename]
            for i, key in enumerate(h[0]):
                if isinstance(h[1][i], str):
                    h[1][i] = h[1][i].strip()
                self.collection[key][row_index] = h[1][i]

    def __getitem__(self, key):
        return self.collection[key]

    def __call__(self, collection):
        ccds = list()
        for filename in collection['filename']:
            ccd = CCDData.read(filename, unit='adu')
            ccds.append(ccd)

        return ccds

    def __str__(self):
        pass

    def __repr__(self):
        pass
