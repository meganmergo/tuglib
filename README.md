# ***tuglib***
Yet another **Astronomical Reduction and Analysis Tool**.

## FitsCollection

```python
from tuglib.reduction import FitsCollection


# Get all fits file on given directory.
images = FitsCollection(location='/Users/oguzhan/tmp/20180901N/')

# Select 'filename' and 'exptime' columns from collection
# where object keyword equal 'BIAS'.
biases = images[images['OBJECT'] == 'BIAS']['filename', 'EXPTIME']
print(biases)
```

```
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
```

```python
# Get ccds object from collection.
ccds = images(biases)

for i, ccd in enumerate(ccds):
    print(biases['filename'][i], '%.3f' % np.average(ccd.data))

```

```
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
```
