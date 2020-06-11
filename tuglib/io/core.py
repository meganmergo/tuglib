#!/usr/bin/env python

__all__ = ['FitsCollection', 'convert_to_ccddata', 'convert_to_fits']


import types
import os
from glob import glob
import fnmatch
from datetime import datetime
from progressbar import ProgressBar

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

    disregard_nan: bool
        If True, any value of nan in the output array will be replaced by
        the readnoise.

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

    stats(include_path=False, **kwargs)
        Creates basic statistical information of the collection.

    Examples
    --------
    >>> from tuglib.io import FitsCollection
    >>>
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

    def __init__(self, location, file_extension='', masks=None, unit=u.adu,
                 gain=None, read_noise=None, disregard_nan=False):

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
        self._disregard_nan = disregard_nan

        self._filenames = list()
        self._filenames_without_path = list()
        self._keywords = None
        self._collection = None

        if self._location:
            self.update()

    def __getitem__(self, key):
        if isinstance(key, np.ndarray):
            return self._collection[key]

        if isinstance(key, str):
            if key != 'filename':
                key = key.upper()

        if isinstance(key, tuple):
            keys = list()
            for item in key:
                if item != 'filename':
                    keys.append(item.upper())
                else:
                    keys.append(item)

            key = keys

        return self._collection[key]

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
        >>> from tuglib.io import FitsCollection
        >>>
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
                    ccd, gain=self._gain, readnoise=self._read_noise,
                    disregard_nan=self._disregard_nan)

                gain_corrected = gain_correct(data_with_deviation, self._gain)

                yield gain_corrected
        else:
            for filename in collection[tmp]['filename']:
                ccd = CCDData.read(filename, unit=self._unit)

                ccd.mask = mask
                ccd = trim_image(ccd, trim)

                yield ccd

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __repr__(self):
        keywords = sorted(self._collection.colnames)
        number_of_keywords = len(keywords)

        number_of_filenames = len(self._filenames)

        objects = list(set(self._collection['OBJECT']))
        number_of_objects = len(objects)

        message = f"Image path: {self._location}\n" \
                  f"{number_of_filenames} images were found.\n\n" \
                  f"Keywords ({number_of_keywords} found):\n" \
                  f"{keywords}\n\n" \
                  f"OBJECT ({number_of_objects} found):\n" \
                  f"{objects}"

        return message

    def __str__(self):
        return 'FITS Image Collection'

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

    def update(self):
        """
        Updates collections.
        """

        if not self._location:
            return

        self._filenames = sorted(
            find_files(self._location, self._file_extension))

        if not self._filenames:
            raise ValueError(f"No images found in '{self._location}'")

        self._filenames_without_path = \
            [os.path.split(filename)[1] for filename in self._filenames]

        self._collection = Table()
        self._collection.add_column(
            Column(data=self._filenames, name='filename'))

        initial_keywords = list(
            set(key for key in fits.getheader(
                self._filenames[0]).keys() if key))

        for keyword in initial_keywords:
            self._collection.add_column(
                Column(data=[None] * len(self._filenames), name=keyword))

        bar = ProgressBar(max_value=len(self._filenames))

        for i, image in enumerate(self._collection['filename']):
            header = fits.getheader(image)
            keywords = list(
                set(key for key in fits.getheader(image).keys() if key))

            for keyword in keywords:
                value = header[keyword]
                if isinstance(value, str):
                    value = value.strip()

                if keyword not in self._collection.colnames:
                    self._collection.add_column(
                        Column(data=[None] * len(self._filenames),
                               name=keyword))

                self._collection[keyword][i] = value

            bar.update(i)

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
            if key == 'filename':
                file_mask = np.array(
                    [fnmatch.fnmatch(filename, kwargs['filename'])
                     for filename in self._filenames_without_path], dtype=bool)

                tmp = tmp & file_mask
            else:
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
            if key == 'filename':
                file_mask = np.array(
                    [fnmatch.fnmatch(filename, kwargs['filename'])
                     for filename in self._filenames_without_path], dtype=bool)

                tmp = tmp & file_mask
            else:
                tmp = tmp & (self._collection[key.upper()] == val)

        filenames = list(self._collection[tmp]['filename'])

        if not filenames:
            raise ValueError('No file found in collection!')

        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        bar = ProgressBar(max_value=len(filenames))

        try:
            ssh.connect(hostname=hostname, username=username,
                        password=password)

            sftp = ssh.open_sftp()

            try:
                sftp.chdir(remote_path)
            except IOError:
                sftp.mkdir(remote_path)
                sftp.chdir(remote_path)

            for i, local_file in enumerate(filenames):
                rel_local_path = os.path.relpath(local_file, self._location)
                remote_file = os.path.join(remote_path, rel_local_path)
                remote_directory = os.path.split(remote_file)[0]

                try:
                    sftp.chdir(remote_directory)
                except IOError:
                    sftp.mkdir(remote_directory)
                    sftp.chdir(remote_directory)

                sftp.put(local_file, remote_file)
                bar.update(i)

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
        >>> from tuglib.io import FitsCollection
        >>>
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
                if key == 'filename':
                    file_mask = np.array(
                        [fnmatch.fnmatch(filename, kwargs['filename'])
                         for filename in self._filenames_without_path],
                        dtype=bool)

                    tmp = tmp & file_mask
                else:
                    tmp = tmp & (self._collection[key.upper()] == val)

        if np.count_nonzero(tmp) == 0:
            yield None

        x = self._collection[tmp]['NAXIS1'][0]
        y = self._collection[tmp]['NAXIS2'][0]
        shape = (y, x)

        mask = None
        if masks is not None:
            mask = make_mask(shape, masks)

        for filename in self._collection[tmp]['filename']:
            ccd = CCDData.read(filename, unit=self._unit,
                               output_verify='silentfix+ignore')

            ccd.mask = mask
            ccd = trim_image(ccd, trim)

            if (self._gain is not None) and (self._read_noise is not None):
                data_with_deviation = create_deviation(
                    ccd, gain=self._gain, readnoise=self._read_noise,
                    disregard_nan=self._disregard_nan)

                gain_corrected = gain_correct(data_with_deviation, self._gain)

                yield gain_corrected
            else:
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

        >>> from tuglib.io import FitsCollection
        >>>
        >>> images = FitsCollection(location='/home/user/data/fits/')
        >>> bias_datas = images.data(OBJECT='BIAS')
        """

        tmp = np.full(len(self._collection), True, dtype=bool)

        if len(kwargs) != 0:
            for key, val in kwargs.items():
                if key == 'filename':
                    file_mask = np.array(
                        [fnmatch.fnmatch(filename, kwargs['filename'])
                         for filename in self._filenames_without_path],
                        dtype=bool)

                    tmp = tmp & file_mask
                else:
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
        >>> from tuglib.io import FitsCollection
        >>>
        >>> images = FitsCollection(location='/home/user/data/fits/')
        >>> bias_headers = images.headers(OBJECT='BIAS')
        """

        tmp = np.full(len(self._collection), True, dtype=bool)

        if len(kwargs) != 0:
            for key, val in kwargs.items():
                if key == 'filename':
                    file_mask = np.array(
                        [fnmatch.fnmatch(filename, kwargs['filename'])
                         for filename in self._filenames_without_path],
                        dtype=bool)

                    tmp = tmp & file_mask
                else:
                    tmp = tmp & (self._collection[key.upper()] == val)

        for filename in self._collection[tmp]['filename']:
            header = fits.getheader(filename)
            yield header

    def stats(self, include_path=False, **kwargs):
        """
        Creates basic statistical information of the collection.

        Parameters
        ----------
        include_path : bool
            Default value is False.

        **kwargs :
            Any additional keywords are used to filter the items returned

        Return
        ------
        'astropy.table.Table'
            Statistical table

        Examples
        --------
        >>> from tuglib.io import FitsCollection
        >>>
        >>> images = FitsCollection(location='/home/user/data/fits/')
        >>> stat = images.stats(OBJECT='BIAS')
        >>> print(stat)
        filename  min    max      mean   median  std     rms
        -------- ----- -------- ------- ------- ------ -------
        bias0001 0.000 2708.000 401.755 400.000 17.024 171.390
        bias0002 0.000 1728.000 398.765 399.000 10.716 167.850
        bias0003 0.000 6627.000 400.096 400.000 11.442 170.973
        bias0004 0.000 5826.000 401.051 401.000 11.348 173.213
        bias0005 0.000 2727.000 401.793 402.000 11.702 174.843
        bias0006 0.000 7634.000 401.558 402.000 11.525 174.370
        bias0007 0.000 3111.000 401.494 402.000 11.929 174.137
        bias0008 0.000 6431.000 401.608 402.000 11.761 174.480
        bias0009 0.000 2658.000 401.361 402.000 10.834 173.927
        bias0010 0.000 3767.000 400.977 401.000 11.020 173.037
        """

        filenames = list()
        mins = list()
        maxs = list()
        means = list()
        medians = list()
        stds = list()
        rmss = list()

        stat = Table()

        tmp = np.full(len(self._collection), True, dtype=bool)

        if len(kwargs) != 0:
            for key, val in kwargs.items():
                if key == 'filename':
                    file_mask = np.array(
                        [fnmatch.fnmatch(filename, kwargs['filename'])
                         for filename in self._filenames_without_path],
                        dtype=bool)

                    tmp = tmp & file_mask
                else:
                    tmp = tmp & (self._collection[key.upper()] == val)

        for filename in self._collection[tmp]['filename']:
            data = fits.getdata(filename)

            if not include_path:
                filename = os.path.split(filename)[1]

            filenames.append(filename)
            mins.append(data.min())
            maxs.append(data.max())
            means.append(data.mean())
            medians.append(np.median(data))
            stds.append(data.std())
            rmss.append(np.sqrt(np.mean(np.square(data.flatten()))))

        stat.add_column(Column(data=filenames, name='filename'))
        stat.add_column(Column(data=mins, name='min'))
        stat.add_column(Column(data=maxs, name='max'))
        stat.add_column(Column(data=means, name='mean'))
        stat.add_column(Column(data=medians, name='median'))
        stat.add_column(Column(data=stds, name='std'))
        stat.add_column(Column(data=rmss, name='rms'))

        stat['min'].info.format = '.3f'
        stat['max'].info.format = '.3f'
        stat['mean'].info.format = '.3f'
        stat['median'].info.format = '.3f'
        stat['std'].info.format = '.3f'
        stat['rms'].info.format = '.3f'

        return stat

    def unique(self, keyword):
        """
        Finds each unique item in keyword.

        Parameters
        ----------

        keyword : str, list or tuple
            The server to connect to.

        Returns
        -------
        list or dict
            Founded unique items.

        Examples
        --------
        >>> from tuglib.io import FitsCollection
        >>>
        >>> images = FitsCollection(location='/home/user/data/fits/')
        >>> print(images.unique('OBJECT'))
        ['Bias', 'Dark', 'Flat', 'Vega']
        """

        if keyword != 'filename':
            keyword = keyword.upper()

        return sorted(list(set(self._collection[keyword])))

    def filter_by_date(self, return_mask=True,
                       date_format='%Y-%m-%dT%H:%M:%S',
                       date_keyword='DATE-OBS', **kwargs):
        """
        Filter collection by date.

        Parameters
        ----------
        return_mask : bool
            If True returns masked array

        date_format : str
            'datetime.datetime' format

        date_keyword : str
            Date keyword to filter. It should be in collection.

        **kwargs :
                Any additional keywords are used to filter the items returned

        Returns
        -------
        'np.array' or 'astropy.table.Table'
            Filtered collection

        Examples
        --------
        >>> from tuglib.io import FitsCollection
        >>>
        >>> images = FitsCollection(location='/home/user/data/fits')
        >>> m = images.filter_by_date(gt='2013-01-10T15:33:00',
                                      lt='2013-01-10T16:22:00')
        >>> print(m)
        array([False, False, False, False, False, False,  True,  True,  True,
                True,  True, False, False, False, False, False, False, False,
               False, False, False, False, False, False, False, False, False,
               False, False]
        """

        if date_keyword not in self._keywords:
            raise ValueError("'date_keyword' error!")

        for key, val in kwargs.items():
            if key not in ['lt', 'gt', 'le', 'ge']:
                raise IOError("It should be used with at least one of the "
                              "'lt', 'gt', 'le' and 'ge' arguments!")

        tmp = np.full(len(self._collection), True, dtype=bool)

        dates = [datetime.strptime(date, date_format)
                 for date in self._collection[date_keyword]]
        t_dates = Column(data=dates, name='dates')

        for key, val in kwargs.items():
            val = datetime.strptime(val, date_format)

            if key == 'lt':
                tmp = tmp & (t_dates < val)
                continue

            if key == 'gt':
                tmp = tmp & (t_dates > val)
                continue

            if key == 'le':
                tmp = tmp & (t_dates <= val)
                continue

            if key == 'ge':
                tmp = tmp & (t_dates >= val)
                continue

            if key.upper() in self._keywords:
                if key == 'filename':
                    file_mask = np.array(
                        [fnmatch.fnmatch(filename, kwargs['filename'])
                         for filename in self._filenames_without_path], dtype=bool)

                    tmp = tmp & file_mask
                else:
                    tmp = tmp & (self._collection[key.upper()] == val)

        if return_mask:
            return tmp

        return self._collection[tmp]
