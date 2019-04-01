#!/usr/bin/env python

__all__ = ['FitsCollection', 'make_mask']

from os import path
from glob import glob

import numpy as np

from astropy.table import Table
from astropy.io import fits
import astropy.units as u

from ccdproc import ImageFileCollection, CCDData, create_deviation,\
    gain_correct, trim_image


# Available file extensions. New extensions can be added in the future.
FILE_EXTENSIONS = ('fit', 'fits', 'fit.gz', 'fits.gz', 'fit.zip', 'fits.zip')


# Another helper method for FitsCollection class.
def get_fits_header(filename):
    """
    Converts 'Header' to list of 'keywords', 'values', 'comments' and 'dtypes'.

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


def make_mask(shape, area):
    """
    Make a mask specific to the given area.

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


class FitsCollection(object):
    """
    FITS Image Collection
    (A ccdproc.ImageFileCollection alternative. It was re written from scratch.)

    It performs recursive fits image search in given directory.

    Parameters
    ----------
    location : str
        Full path directory that include images.
        Example: '/home/user/data/20190325'

    file_extension : None, str, tuple or list
        Image extensions. Should be one of 'fit', 'fits', 'fit.gz',
        'fits.gz', 'fit.zip', 'fits.zip'

    gain : None or float
        'Gain' value. Default type is 'None'.

    read_noise : None or float
        'Read Noise' value. Default type is 'None'.

    unit : astropy.units
        'Unit' value. Default type is 'astropy.units.adu'.

    Methods
    -------
    ccds(**kwargs)
        Return 'Data' objects from collection. If any positional
        arguments like 'fits' header keywords, it filters result.

    datas(**kwargs)
        Return 'np.array' (fits.data) objects from collection. If any positional
        arguments like 'fits' header keywords, it filters result.

    headers(**kwargs)
        Return 'fits.header' objects from collection. If any positional
        arguments like 'fits' header keywords, it filters result.

    Examples
    --------

    >>> # Get all fits file from directory.
    >>> images = FitsCollection(
            location='/Users/oguzhan/tmp/20180901N/',
            gain=0.57, read_noise=4.11)
    >>>
    >>> # Get 'ccd' objects from collection where 'object' keyword equal 'BIAS'.
    >>> biases = images.ccds(OBJECT='BIAS')
    >>>
    >>> # Filtered search example.
    >>> # Select 'filename' and 'exptime' columns from collection
    >>> # where 'object' keyword equals 'FLAT' and
    >>> # 'filter' keyword equals 'W1:03 V W2:00 Empty' and
    >>> # 'exptime' keyword less than 0.06 seconds.
    >>> f1 = images['OBJECT'] == 'FLAT'
    >>> f2 = images['FILTER'] == 'W1:03 V W2:00 Empty'
    >>> f3 = images['EXPTIME'] < 0.06
    >>>
    >>> filtered_flats = images[f1 & f2 & f3]['filename', 'EXPTIME']
    >>> print(filtered_flats)
                         filename                     EXPTIME
                          str64                       float64
    ------------------------------------------------- -------
    /Users/oguzhan/tmp/20180901N/BDF/FLAT_0018_V.fits    0.05
    /Users/oguzhan/tmp/20180901N/BDF/FLAT_0019_V.fits    0.05
    """

    def __init__(self, location, file_extension='fits', masks=None,
                 gain=None, read_noise=None, unit=u.adu):

        if not isinstance(location, str):
            raise TypeError(
                "'location' should be a 'str' object.")

        if not isinstance(file_extension, (str, tuple, list)):
            raise TypeError(
                "'file_extension' should be "
                "'None', 'str', 'tuple' or 'list' object.")

        if file_extension not in FILE_EXTENSIONS:
            raise ValueError(
                "'file_extension' should be "
                "'fit, fits, fit.gz, fits.gz, fit.zip, fits.zip' type.")

        if not isinstance(gain, (type(None), float)):
            raise TypeError(
                "'gain' should be 'None' or 'float' object.")

        if not isinstance(read_noise, (type(None), float)):
            raise TypeError(
                "'read_noise' should be 'None' or 'float' object.")

        self._location = location
        self._file_extension = file_extension

        self._masks = masks

        self._gain = gain
        if gain is not None:
            self._gain = gain * u.electron / u.adu

        self._read_noise = read_noise
        if read_noise is not None:
            self._read_noise = read_noise * u.electron

        self._unit = unit

        self._filenames = list()
        self._keywords = None
        self._collection = None

        self._prepare()

    @property
    def location(self):
        return self._location

    @property
    def gain(self):
        return self._gain

    @property
    def read_noise(self):
        return self._read_noise

    @property
    def unit(self):
        return self._unit

    @property
    def filenames(self):
        return self._filenames

    @property
    def keywords(self):
        return self._keywords

    @property
    def collection(self):
        return self._collection

    def _prepare(self):
        """
        Only for internal use.
        """
        directory = list()

        if self._file_extension is None:
            for ext in FILE_EXTENSIONS:
                directory += path.join(self._location, '**', '*.' + ext)
        elif isinstance(self._file_extension, (tuple, list)):
            for ext in self._file_extension:
                directory += path.join(self._location, '**', '*.' + ext)
        else:  # str
            directory = path.join(
                self._location, '**', '*.' + self._file_extension)

        self._filenames = sorted(glob(directory, recursive=True))

        db = dict()
        for filename in self._filenames:
            h = get_fits_header(filename)
            db[filename] = h

        cs = dict()
        for filename in self._filenames:
            for i, key in enumerate(db[filename][0]):
                cs[db[filename][0][i]] = db[filename][3][i]

        names = list(cs.keys())
        names = ['filename'] + names
        dtypes = list(cs.values())
        dtypes = ['U64'] + dtypes

        self._collection = Table(names=names, dtype=dtypes)

        for filename in self._filenames:
            self._collection.add_row()
            row_index = len(self._collection) - 1
            self._collection['filename'][row_index] = filename

            h = db[filename]
            for i, key in enumerate(h[0]):
                if isinstance(h[1][i], str):
                    h[1][i] = h[1][i].strip()
                self._collection[key][row_index] = h[1][i]

        self._keywords = self._collection.colnames

    def ccds(self, masks=None, trim=None, **kwargs):
        """
        Get 'CCDData' objects from collection.

        Parameters
        ----------
        masks : list, optional
            A list of string formated image mask section.

        trim : str, optional
            Trim section.

        **kwargs :
            Any additional keywords are used to filter the items returned.

        Returns
        -------
        ccds: list of CCDData

        Examples
        --------

        >>> images = FitsCollection(location='/home/user/data/fits/')
        >>> biases = images.ccds(OBJECT='BIAS')
        """

        if masks is not None:
            if not isinstance(masks, list):
                raise TypeError("'masks' should be a 'list' object.")

        if trim is not None:
            if not isinstance(trim, str):
                raise TypeError("'trim' should be a 'str' object.")

        records = list()

        tmp = np.full(len(self._collection), True, dtype=bool)

        if len(kwargs) != 0:
            for key, val in kwargs.items():
                tmp = tmp & (self._collection[key] == val)

        x = self._collection[tmp]['NAXIS1'][0]
        y = self._collection[tmp]['NAXIS2'][0]
        shape = (y, x)

        mask = None
        if masks is not None:
            mask = make_mask(shape, masks[0])
            for m in masks[1:]:
                tmp_mask = make_mask(shape, m)
                mask |= tmp_mask

        if (self._gain is not None) and (self._read_noise is not None):
            for filename in self._collection[tmp]['filename']:
                ccd = CCDData.read(filename, unit=self._unit)

                # if mask is not None:
                ccd.mask = mask
                ccd = trim_image(ccd, trim)

                data_with_deviation = create_deviation(
                    ccd, gain=self._gain, readnoise=self._read_noise)

                gain_corrected = gain_correct(data_with_deviation, self._gain)

                records.append(gain_corrected)

            return records

        for filename in self._collection[tmp]['filename']:
            ccd = CCDData.read(filename, unit=self._unit)
            records.append(ccd)

        return records

    def datas(self, **kwargs):
        """
        Get 'numpy.ndarray' objects from collection.

        Parameters
        ----------
        **kwargs :
            Any additional keywords are used to filter the items returned

        Returns
        -------
        datas: list of 'numpy.ndarray'

        Examples
        --------

        >>> images = FitsCollection(location='/home/user/data/fits/')
        >>> bias_datas = images.datas(OBJECT='BIAS')
        """
        records = list()

        tmp = np.full(len(self._collection), True, dtype=bool)

        if len(kwargs) != 0:
            for key, val in kwargs.items():
                tmp = tmp & (self._collection[key] == val)

        for filename in self._collection[tmp]['filename']:
            data = fits.getdata(filename)
            records.append(data)

        return records

    def headers(self, **kwargs):
        """
        Get 'astropy.io.fits.header.Header' objects from collection.

        Parameters
        ----------
        **kwargs :
            Any additional keywords are used to filter the items returned

        Returns
        -------
        datas: list of 'astropy.io.fits.header.Header'

        Examples
        --------

        >>> images = FitsCollection(location='/home/user/data/fits/')
        >>> bias_headers = images.headers(OBJECT='BIAS')
        """
        records = list()

        tmp = np.full(len(self._collection), True, dtype=bool)

        if len(kwargs) != 0:
            for key, val in kwargs.items():
                tmp = tmp & (self._collection[key] == val)

        for filename in self._collection[tmp]['filename']:
            data = fits.getheader(filename)
            records.append(data)

        return records

    def _hdus(self, **kwargs):
        records = list()

        tmp = np.full(len(self._collection), True, dtype=bool)

        if len(kwargs) != 0:
            for key, val in kwargs.items():
                tmp = tmp & (self._collection[key] == val)

        for filename in self._collection[tmp]['filename']:
            hdu = fits.open(filename)
            records.append(hdu)

        return records

    def __getitem__(self, key):
        return self._collection[key]

    def __call__(self, collection=None, **kwargs):
        if collection is not None:
            ccds = list()

            for filename in collection['filename']:
                ccd = CCDData.read(filename, unit=self._unit)
                ccds.append(ccd)

            return ccds

        return self.ccds(**kwargs)


# It doesn't work well. Just for proof of concept.
# Helper class for ccdproc.ImageFileCollection
class FitsCollectionTest(ImageFileCollection):

    def __init__(self, **kwargs):
        if 'location' in kwargs:
            directory = path.join(kwargs['location'], '**', '*.*')

            images = glob(directory, recursive=True)
            # kwargs['filenames'] = images
            kwargs['filenames'] = \
                [image.split(kwargs['location'])[1]
                 for image in images]
            print(kwargs['filenames'])

        super(FitsCollectionTest, self).__init__(**kwargs)
