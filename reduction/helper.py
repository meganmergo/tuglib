#!/usr/bin/env python

__all__ = ['FitsCollectionTest']

from os import path
from glob import glob

from astropy import units as u
from astropy.table import Table
from astropy.io import fits

from ccdproc import ImageFileCollection, CCDData


FILE_EXTENSIONS = ('fit', 'fits', 'fit.gz', 'fits.gz', 'fit.zip', 'fits.zip')


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
    """FITS Image Collection.
    . . .
    """

    def __init__(self, location=None, file_extension='fits', **kwargs):

        if not isinstance(location, (None, str)):
            raise TypeError(
                "'location' should be 'None' or 'str' object.")

        if not isinstance(file_extension, (str, tuple, list)):
            raise TypeError(
                "'file_extension' should be 'str', 'tuple' or 'list' object.")

        if file_extension not in FILE_EXTENSIONS:
            raise ValueError(
                "'file_extension' should be "
                "'fit, fits, fit.gz, fits.gz, fit.zip, fits.zip'.")

        self.location = location
        self._file_extension = file_extension
        self.filenames = list()
        self.collection = Table()

        if location is not None:
            self.init()

    def init(self):
        directory = path.join(self.location, '**', '*.*')
        self.filenames = glob(directory, recursive=True)

        for filename in self.filenames:
            header = fits.getheader(filename)
            # header = ccd.meta.cards


    def __str__(self):
        pass

    def __repr__(self):
        pass
