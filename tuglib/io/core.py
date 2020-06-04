#!/usr/bin/env python

__all__ = ['FitsCollection', 'convert_to_ccddata', 'convert_to_fits']


import types
import os
from glob import glob
from progress.bar import Bar

import paramiko

import numpy as np

from astropy.table import Table, vstack, Column
from astropy.io import fits
import astropy.units as u

from ccdproc import CCDData, create_deviation, gain_correct, trim_image

from ..reduction.helper import make_mask
from .helper import find_files

# Available file extensions. New extensions can be added in the future.
FILE_EXTENSIONS = ('fit', 'fits', 'fit.gz', 'fits.gz', 'fit.zip', 'fits.zip',
                   'fts', 'fts.gz', 'fts.zip', '')


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
            os.mkdir(location)
        except FileExistsError:
            pass
    else:
        location = ''

    if isinstance(images, CCDData):
        if not isinstance(filenames, str):
            raise TypeError("'filenames' should be a 'str' object.")

        output = filenames.split(os.sep)[-1]

        if prefix is not None:
            output = prefix + output

        if suffix is not None:
            tmp = output.split('.')

            output = tmp[0] + suffix

            if len(tmp) > 1:
                output = output + '.' + tmp[-1]

        output = os.path.join(location, output)
        images.write(output, overwrite=True, output_verify='silentfix+ignore')
    else:
        for i, ccd in enumerate(images):
            output = filenames[i].split(os.sep)[-1]

            if prefix is not None:
                output = prefix + output

            if suffix is not None:
                tmp = output.split('.')

                output = tmp[0] + suffix

                if len(tmp) > 1:
                    output = output + '.' + tmp[-1]

            output = os.path.join(location, output)
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

    upload(hostname, username, password, remote_path, **kwargs)
        Upload 'fits' file to remote computer with SFTP.

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

    def __init__(self, location, file_extension='', masks=None,
                 unit=u.adu, gain=None, read_noise=None):

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

        if self._location:
            self._prepare()

    def __getitem__(self, key):
        return self._collection[key.upper()]

    def __add__(self, collection):
        if not (self._keywords == collection.keywords):
            raise ValueError(
                'The header keywords of the two collections must be the same.')

        if self._gain != collection.gain:
            raise ValueError(
                'The gain of the two collections must be the same.')

        if self._read_noise != collection.read_noise:
            raise ValueError(
                'The read_noise of the two collections must be the same.')

        if self._unit != collection.unit:
            raise ValueError(
                'The unit of the two collections must be the same.')

        c = FitsCollection('')

        c._concat(self._collection, collection.collection.copy(),
                  self._location, collection.location,
                  self._filenames, collection.filenames.copy(),
                  self._keywords, self._gain, self._read_noise, self._unit)

        return c

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
    def collection(self):
        return self._collection

    def _prepare(self):
        """
        Only for internal use.
        """
        self._filenames = sorted(
            find_files(self._location, self._file_extension))

        self._collection = Table()
        self._collection.add_column(
            Column(data=self._filenames, name='filename'))

        initial_keywords = list(
            set(key for key in fits.getheader(
                self._filenames[0]).keys() if key))

        for keyword in initial_keywords:
            self._collection.add_column(
                Column(data=[None] * len(self._filenames), name=keyword))

        bar = Bar('Creating FitsCollection', max=len(self._filenames),
                  fill='#', suffix='%(percent).1f%% - %(eta)ds')

        for i, image in enumerate(self._collection['filename']):
            header = fits.getheader(image)
            keywords = list(
                set(key for key in fits.getheader(image).keys() if key))

            for keyword in keywords:
                try:
                    self._collection[keyword][i] = header[keyword]
                except KeyError:
                    self._collection.add_column(
                        Column(data=[None] * len(self._filenames),
                               name=keyword))

                    self._collection[keyword][i] = header[keyword]

            bar.next()

        bar.finish()

        self._keywords = self._collection.colnames

    def _concat(self, collection_1, collection_2, location_1, location_2,
                filenames_1, filenames_2, keywords_1, gain, read_noise, unit):
        """
        Only for internal use.
        """

        c = [collection_1, collection_2]
        collection = vstack(c)

        filenames = filenames_1 + filenames_2

        self._collection = collection

        self._location = [location_1, location_2]
        self._filenames = filenames
        self._keywords = keywords_1.copy()

        self._gain = gain
        self._read_noise = read_noise
        self._unit = unit

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
            tmp = tmp & (self._collection[key.upper()] == val)

        if np.count_nonzero(tmp) == 0:
            return list()

        files = list(self._collection[tmp]['filename'])

        if include_path:
            return files

        files = [file.split(os.sep)[-1] for file in files]

        return files

    def upload(self, hostname, username, password, remote_path, **kwargs):
        """
        Upload 'fits' file to remote computer with SFTP.

        Parameters
        ----------

        hostname : str
            The server to connect to.

        username : str
            The username to authenticate as.

        password : str
            Used for password authentication.

        remote_path : str
            The destination path on the SFTP server.

        Returns
        -------

        Bool

        Examples
        --------

        >>> from tuglib.io import FitsCollection
        >>>
        >>> c = FitsCollection('/home/user/images')
        >>> c.upload(hostname='192.168.0.1', username='guest',
                     password='1234', remote_path='/home/guest/works')
        """

        if not isinstance(hostname, str):
            raise TypeError("'hostname' should be a str object.")

        if not isinstance(username, str):
            raise TypeError("' should be a str object.")

        if not isinstance(password, str):
            raise TypeError("'password' should be a str object.")

        if not isinstance(remote_path, str):
            raise TypeError("'remote_path' should be a str object.")

        tmp = np.full(len(self._collection), True, dtype=bool)

        for key, val in kwargs.items():
            tmp = tmp & (self._collection[key.upper()] == val)

        filenames = list(self._collection[tmp]['filename'])

        if not filenames:
            raise ValueError('No file found in collection!')

        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        bar = Bar('Uploading', max=len(filenames),
                  fill='#', suffix='%(percent).1f%% - %(eta)ds')

        try:
            ssh.connect(hostname=hostname, username=username,
                        password=password)

            sftp = ssh.open_sftp()

            try:
                sftp.chdir(remote_path)
            except IOError:
                sftp.mkdir(remote_path)
                sftp.chdir(remote_path)

            for local_file in filenames:
                rel_local_path = os.path.relpath(local_file, self._location)
                remote_file = os.path.join(remote_path, rel_local_path)
                remote_directory = os.path.split(remote_file)[0]

                try:
                    sftp.chdir(remote_directory)
                except IOError:
                    sftp.mkdir(remote_directory)
                    sftp.chdir(remote_directory)

                sftp.put(local_file, remote_file)
                bar.next()

            sftp.close()
            ssh.close()
        except IOError as e:
            raise IOError(e)

        bar.finish()

        return True

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
                tmp = tmp & (self._collection[key.upper()] == val)

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
                tmp = tmp & (self._collection[key.upper()] == val)

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
                tmp = tmp & (self._collection[key.upper()] == val)

        for filename in self._collection[tmp]['filename']:
            header = fits.getheader(filename)
            yield header

