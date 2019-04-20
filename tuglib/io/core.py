#!/usr/bin/env python

__all__ = ['FitsCollection', 'convert_to_ccddata', 'convert_to_fits']


import types
from os import path, sep, mkdir
from glob import glob

import numpy as np

from astropy.table import Table
from astropy.io import fits
import astropy.units as u

from ccdproc import CCDData, create_deviation, gain_correct, trim_image

from ..reduction.helper import make_mask

# Available file extensions. New extensions can be added in the future.
FILE_EXTENSIONS = ('fit', 'fits', 'fit.gz', 'fits.gz', 'fit.zip', 'fits.zip',
                   'fts', 'fts.gz', 'fts.zip')


def convert_to_ccddata(images, gain=None, read_noise=None):
    """
    Convert 'fits' file to 'ccdproc.CCCData' object.

    Parameters
    ----------
    images : str or list of str.
        Images to converted.

    gain : float
        Gain (u.electron / u.adu)

    read_noise : float
        Read Noise (u.electron).

    Yields
    ------
    'ccdproc.CCDData'
        yield the next 'ccdproc.CCDData'.

    Examples
    --------
    >>> from tuglib.io import convert_to_ccddata
    >>> from glob import glob
    >>>
    >>> images = glob('/home/user/data/image*.fits')
    >>>
    >>> ccds = convert_to_ccddata(images, gain=0.37, read_noise=4.11)
    """

    if not isinstance(images, (str, list)):
        raise TypeError(
            "'images' should be 'str' or 'list' object.")

    if not isinstance(gain, (type(None), float)):
        raise TypeError("'gain' should be 'None' or 'float' object.")

    if not isinstance(read_noise, (type(None), float)):
        raise TypeError("'read_noise' should be 'None' or 'float' object.")

    if isinstance(images, str):
        images = [images]

    unit = u.adu

    if gain is not None:
        gain = gain * u.electron / u.adu

    if read_noise is not None:
        read_noise = read_noise * u.electron

    for image in images:
        ccd = CCDData.read(image, unit=unit, output_verify='silentfix+ignore')

        if (gain is not None) and (read_noise is not None):
            data_with_deviation = create_deviation(
                ccd, gain=gain, readnoise=read_noise)

            gain_corrected = gain_correct(data_with_deviation, gain)

            yield gain_corrected
        else:
            yield ccd


def convert_to_fits(images, filenames, location=None,
                    prefix=None, suffix=None):
    """
    Convert 'ccdproc.CCCData' to 'fits' file.

    Parameters
    ----------
    images : 'CCDData', list of 'CCDData' or generator.
        Images to converted.

    filenames : str or list of str
        Output filenames.

    location : optional, str
        Path where to save the 'fits' files. Default path is current path.

    prefix : optional, str
        Output file prefix.

    suffix : str
        Output file suffix.

    Returns
    -------
    'fits' file.

    Examples
    --------
    >>> from tuglib.io import FitsCollection, convert_to_fits
    >>>
    >>> c = FitsCollection('/home/user/data')
    >>>
    >>> bias_ccds = c.ccds(OBJECT='BIAS')
    >>> filenames = c[c['OJBECT'] == 'BIAS']['filename']
    >>> convert_to_fits(bias_ccds, filenames=filenames, prefix='test_')
    """

    if not isinstance(images, (CCDData, list, types.GeneratorType)):
        raise TypeError(
            "'images' should be 'CCDData', 'list' or 'generator'.")

    if not isinstance(filenames, (str, list)):
        raise TypeError("'filenames' should be 'str' or 'list' object.")

    if not isinstance(location, (type(None), str)):
        raise TypeError("'prefix' should be a 'str' object.")

    if not isinstance(prefix, (type(None), str)):
        raise TypeError("'prefix' should be a 'str' object.")

    if not isinstance(suffix, (type(None), str)):
        raise TypeError("'suffix' should be a 'int' object.")

    if location is not None:
        try:
            mkdir(location)
        except FileExistsError:
            pass
    else:
        location = ''

    if isinstance(images, CCDData):
        if not isinstance(filenames, str):
            raise TypeError("'filenames' should be a 'str' object.")

        output = filenames.split(sep)[-1]

        if prefix is not None:
            output = prefix + output

        if suffix is not None:
            tmp = output.split('.')

            output = tmp[0] + suffix

            if len(tmp) > 1:
                output = output + '.' + tmp[-1]

        output = path.join(location, output)
        images.write(output, overwrite=True, output_verify='silentfix+ignore')
    else:
        for i, ccd in enumerate(images):
            output = filenames[i].split(sep)[-1]

            if prefix is not None:
                output = prefix + output

            if suffix is not None:
                tmp = output.split('.')

                output = tmp[0] + suffix

                if len(tmp) > 1:
                    output = output + '.' + tmp[-1]

            output = path.join(location, output)
            ccd.write(output, overwrite=True, output_verify='silentfix+ignore')


# Helper method for FitsCollection class.
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
        Return 'ccdproc.CCDData' objects from collection. If any positional
        arguments like 'fits' header keywords, it filters result.

    data(**kwargs)
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
        dtypes = ['U256'] + dtypes

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

    def files_filtered(self, include_path=False, **kwargs):
        """
        Determine files whose keywords have listed values.

        Parameters
        ----------
        include_path : bool
            If is True, returned files include full path.
            Default is False.

        **kwargs :
            Any additional keywords are used to filter the items returned.

        Returns
        -------
        list of str
            Filtered file names from collection.

        Examples
        --------

        >>> from tuglib.io import FitsCollection
        >>>
        >>> c = FitsCollection('/home/user/data')
        >>> files = c.files_filtered(OBJECT='BIAS')
        """

        tmp = np.full(len(self._collection), True, dtype=bool)

        for key, val in kwargs.items():
            tmp = tmp & (self._collection[key] == val)

        if np.count_nonzero(tmp) == 0:
            return list()

        files = list(self._collection[tmp]['filename'])

        if include_path:
            return files

        files = [file.split(sep)[-1] for file in files]

        return files

    def ccds(self, masks=None, trim=None, **kwargs):
        """
        Generator that yields each 'ccdproc.CCDData' objects in the collection.

        Parameters
        ----------
        masks : str, list of str or optional
            Area to be masked.

        trim : str or optional
            Trim section.

        **kwargs :
            Any additional keywords are used to filter the items returned.

        Yields
        ------
        'ccdproc.CCDData'
            yield the next 'ccdproc.CCDData' in the collection.

        Examples
        --------

        >>> mask = '[:, 1000:1046]'
        >>> trim = '[100:1988, :]'
        >>> images = FitsCollection(
                location='/home/user/data/fits/', gain=0.57, read_noise=4.11)
        >>> biases = images.ccds(OBJECT='BIAS', masks=mask, trim=trim)
        """

        if masks is not None:
            if not isinstance(masks, (str, list, type(None))):
                raise TypeError(
                    "'masks' should be 'str', 'list' or 'None' object.")

        if trim is not None:
            if not isinstance(trim, str):
                raise TypeError("'trim' should be a 'str' object.")

        tmp = np.full(len(self._collection), True, dtype=bool)

        if len(kwargs) != 0:
            for key, val in kwargs.items():
                tmp = tmp & (self._collection[key] == val)

        if np.count_nonzero(tmp) == 0:
            yield None

        x = self._collection[tmp]['NAXIS1'][0]
        y = self._collection[tmp]['NAXIS2'][0]
        shape = (y, x)

        mask = None
        if masks is not None:
            mask = make_mask(shape, masks)

        if (self._gain is not None) and (self._read_noise is not None):
            for filename in self._collection[tmp]['filename']:
                ccd = CCDData.read(filename, unit=self._unit,
                                   output_verify='silentfix+ignore')

                ccd.mask = mask
                ccd = trim_image(ccd, trim)

                data_with_deviation = create_deviation(
                    ccd, gain=self._gain, readnoise=self._read_noise)

                gain_corrected = gain_correct(data_with_deviation, self._gain)

                yield gain_corrected
        else:
            for filename in self._collection[tmp]['filename']:
                ccd = CCDData.read(filename, unit=self._unit,
                                   output_verify='silentfix+ignore')

                ccd.mask = mask
                ccd = trim_image(ccd, trim)

                yield ccd

    def data(self, **kwargs):
        """
        Generator that yields each 'numpy.ndarray' objects in the collection.

        Parameters
        ----------
        **kwargs :
            Any additional keywords are used to filter the items returned

        Yields
        ------
        'numpy.ndarray'
            yield the next 'numpy.ndarray' in the collection.

        Examples
        --------

        >>> images = FitsCollection(location='/home/user/data/fits/')
        >>> bias_datas = images.data(OBJECT='BIAS')
        """

        tmp = np.full(len(self._collection), True, dtype=bool)

        if len(kwargs) != 0:
            for key, val in kwargs.items():
                tmp = tmp & (self._collection[key] == val)

        for filename in self._collection[tmp]['filename']:
            yield fits.getdata(filename)

    def headers(self, **kwargs):
        """
        Generator that yields each 'astropy.io.fits.header.Header'
        objects in the collection.

        Parameters
        ----------
        **kwargs :
            Any additional keywords are used to filter the items returned

        Yields
        ------
        'astropy.io.fits.header.Header'
            yield the next 'astropy.io.fits.header.Header' in the collection.

        Examples
        --------

        >>> images = FitsCollection(location='/home/user/data/fits/')
        >>> bias_headers = images.headers(OBJECT='BIAS')
        """

        tmp = np.full(len(self._collection), True, dtype=bool)

        if len(kwargs) != 0:
            for key, val in kwargs.items():
                tmp = tmp & (self._collection[key] == val)

        for filename in self._collection[tmp]['filename']:
            header = fits.getheader(filename)
            yield header

    def __getitem__(self, key):
        return self._collection[key]

    def __call__(self, collection=None, masks=None, trim=None, **kwargs):
        """
        Generator that yields each 'ccdproc.CCDData' objects in the collection.

        Parameters
        ----------
        collection : 'FitsCollection.collection' or optional
            Filtered collection.

        masks : str, list of str or optional
            Area to be masked.

        trim : str or optional
            Trim section.

        **kwargs :
            Any additional keywords are used to filter the items returned.

        Yields
        ------
        'ccdproc.CCDData'
            yield the next 'ccdproc.CCDData' in the collection.

        Examples
        --------

        >>> mask = '[:, 1000:1046]'
        >>> trim = '[100:1988, :]'
        >>>
        >>> images = FitsCollection(
                location='/home/user/data/fits/', gain=0.57, read_noise=4.11)
        >>>
        >>> query = images['EXPTIME'] == 100.0
        >>> sub_collections = images[query]
        >>>
        >>> ccds = images(sub_collections, masks=mask, trim=trim)
        """

        if masks is not None:
            if not isinstance(masks, (str, list, type(None))):
                raise TypeError(
                    "'masks' should be 'str', 'list' or 'None' object.")

        if trim is not None:
            if not isinstance(trim, str):
                raise TypeError("'trim' should be a 'str' object.")

        if collection is None:
            return self.ccds(masks=masks, trim=trim, **kwargs)

        tmp = np.full(len(collection), True, dtype=bool)

        if len(kwargs) != 0:
            for key, val in kwargs.items():
                tmp = tmp & (collection[key] == val)

        x = collection[tmp]['NAXIS1'][0]
        y = collection[tmp]['NAXIS2'][0]
        shape = (y, x)

        mask = None
        if masks is not None:
            mask = make_mask(shape, masks)

        if (self._gain is not None) and (self._read_noise is not None):
            for filename in collection[tmp]['filename']:
                ccd = CCDData.read(filename, unit=self._unit)

                ccd.mask = mask
                ccd = trim_image(ccd, trim)

                data_with_deviation = create_deviation(
                    ccd, gain=self._gain, readnoise=self._read_noise)

                gain_corrected = gain_correct(data_with_deviation, self._gain)

                yield gain_corrected
        else:
            for filename in collection[tmp]['filename']:
                ccd = CCDData.read(filename, unit=self._unit)

                ccd.mask = mask
                ccd = trim_image(ccd, trim)

                yield ccd