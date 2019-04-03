# ***tuglib***
Yet another **Astronomical Reduction and Analysis Tool**.

## FitsCollection

```python
from tuglib.reduction import FitsCollection


# Get all fits file from directory.
images = FitsCollection(
    location='/Users/oguzhan/tmp/20180901N/', gain=0.57, read_noise=4.11)

# Get 'ccd' objects from collection where 'object' keyword equal 'BIAS'.
biases = images.ccds(OBJECT='BIAS')

# Filtered search example.
# Select 'filename' and 'exptime' columns from collection
# where 'object' keyword equals 'FLAT' and
# 'filter' keyword equals 'W1:03 V W2:00 Empty' and
# 'exptime' keyword less than 0.06 seconds.
f1 = images['OBJECT'] == 'FLAT'
f2 = images['FILTER'] == 'W1:03 V W2:00 Empty'
f3 = images['EXPTIME'] < 0.06

filtered_flats = images[f1 & f2 & f3]['filename', 'EXPTIME']
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

## Image Combine

#### Generic Combine

```python
from tuglib.reduction import FitsCollection, image_combine


images = FitsCollection('/home/user/data/m31')

master_image = image_combine(images, method='sum', output='master.fits')
```

#### Bias Combine

```python
from tuglib.reduction import FitsCollection, bias_combine


path = '/home/user/data'
masks = ['[:, 1023:1025]', '[:1023, 56:58]']
trim = '[:, 24:2023]'

images = FitsCollection(location=path, gain=0.57, read_noise=4.11)

master_bias = bias_combine(images, method='median', output='master_bias.fits',
                           masks=masks, trim=trim, OBJECT='BIAS')
```
