#!/usr/bin/env python

__all__ = ['FitsCollection']

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
# Helper class for ccdproc.ImageFileCollection
class FitsCollectionTest(ImageFileCollection):

    def __init__(self, **kwargs):
        if 'location' in kwargs:
            directory = path.join(kwargs['location'], '**', '*.*')

            filenames = glob(directory, recursive=True)
            kwargs['filenames'] = filenames

        super(FitsCollectionTest, self).__init__(**kwargs)


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

    Examples
    --------

    >>> # Get all fits file on given directory.
    >>> images = FitsCollection(location='/Users/oguzhan/tmp/20180901N/')
    >>>
    >>> # Select 'filename' and 'exptime' columns from collection
    >>> # where object keyword equal 'BIAS'.
    >>> biases = images[images['OBJECT'] == 'BIAS']['filename', 'EXPTIME']
    >>> print(biases)
                        filename                     EXPTIME
    ------------------------------------------------ -------
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0001.fits     0.0
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0002.fits     0.0
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0003.fits     0.0
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0004.fits     0.0
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0005.fits     0.0
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0006.fits     0.0
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0007.fits     0.0
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0008.fits     0.0
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0009.fits     0.0
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0010.fits     0.0
    >>>
    >>> # Get ccds object from collection.
    >>> ccds = images(biases)
    >>>
    >>> for i, ccd in enumerate(ccds):
    >>>     print(biases['filename'][i], '%.3f' % np.average(ccd.data))
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0001.fits 496.991
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0002.fits 496.999
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0003.fits 497.006
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0004.fits 496.996
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0005.fits 497.001
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0006.fits 497.004
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0007.fits 497.000
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0008.fits 496.994
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0009.fits 496.987
    /Users/oguzhan/tmp/20180901N/BDF/Bias__0010.fits 496.980
    """

    def __init__(self, location, file_extension=None):

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
