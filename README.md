# ***tuglib***
Yet another **Astronomical Reduction and Analysis Tool**.

## FitsCollection

### Basic usage

```python
from tuglib.io import FitsCollection


# Get all fits file from directory.
images = FitsCollection(
    location='/Users/oguzhan/tmp/20180901N/', gain=0.57, read_noise=4.11)

# Get 'ccd' objects from collection where 'object' keyword equal 'BIAS'.
biases = images.ccds(OBJECT='BIAS')

# Filtered search example.
# Select 'filename' and 'exptime' columns from collection
# where 'object' or 'imagetyp' keywords equals 'FLAT' and
# 'filter' keyword equals 'W1:03 V W2:00 Empty' and
# 'exptime' keyword less than 0.06 seconds.

f0 = images['OBJECT'] == 'FLAT'
f1 = images['IMAGETYP'] == 'FLAT
f2 = images['FILTER'] == 'W1:03 V W2:00 Empty'
f3 = images['EXPTIME'] < 0.06

filtered_flats = images[(f0 | f1) & f2 & f3]['filename', 'EXPTIME']
print(filtered_flats)
```

```
                     filename                     EXPTIME
                      str64                       float64
------------------------------------------------- -------
/Users/oguzhan/tmp/20180901N/BDF/FLAT_0018_V.fits    0.05
/Users/oguzhan/tmp/20180901N/BDF/FLAT_0019_V.fits    0.05
```

```python
# Get 'ccd' objects from 'filtered_flats' collection.
flats = images(filtered_flats)
```

#### Uploading 'fits' files to remote computer

```python
from tuglib.io import FitsCollection

images = FitsCollection(location='/home/user/datas/20190621')

images.upload(hostname='192.168.0.1', username='guest', password='1234',
              remote_path='/home/guest/works', OBJECT='BIAS')
```

## Image Combine

#### Generic Combine

```python
from tuglib.io import FitsCollection
from tuglib.reduction import image_combine


images = FitsCollection('/home/user/data/m31')

ccds = images.ccds(OBJECT='M31_R')
master_image = image_combine(ccds, method='sum', output='master.fits')
```

#### Bias, Dark and Flat Combine

```python
from tuglib.io import FitsCollection
from tuglib.reduction import bias_combine, dark_combine, flat_combine


path = '/home/user/data/20181026'
masks = ['[:, 1023:1025]', '[:1023, 56:58]']
trim = '[:, 24:2023]'

images = FitsCollection(location=path, gain=0.57, read_noise=4.11)

bias_ccds = images.ccds(OBJECT='BIAS', trim=trim, masks=masks)
dark_ccds = images.ccds(OBJECT='DARK', EXPTIME=10, trim=trim, masks=masks)
flat_ccds = images.ccds(OBJECT='FLAT', FILTER='V', trim=trim, masks=masks)

master_bias = bias_combine(bias_ccds, method='median')
master_dark = dark_combine(dark_ccds, master_bias, method='median')
master_flat = flat_combine(flat_ccds, master_bias, master_dark, method='median')
```

#### Basic Reduction (Bias, Dark and Flat Corrections)

```python
from tuglib.reduction import ccdproc


sci_ccds = images.ccds(OBJECT='Star', FILTER='V', trim=trim, masks=masks)

# Yield a generator which point reduced images.
reduced_ccds = ccdproc(sci_ccds, master_bias, master_dark, master_flat)
```
