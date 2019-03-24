# tuglib
Yet another astronomical library 

##FitsCollection

```python
from tuglib.reduction import FitsCollection


c = FitsCollection(location='/Users/oguzhan/tmp/20180901N/BDF',
                   file_extension='fits')

print(c.collection[c.collection['EXPTIME'] > 0.1]['filename', 'EXPTIME'])

```

```
                     filename                     EXPTIME
------------------------------------------------- -------
/Users/oguzhan/tmp/20180901N/BDF/FLAT_0001_U.fits     0.3
/Users/oguzhan/tmp/20180901N/BDF/FLAT_0002_U.fits     0.2
/Users/oguzhan/tmp/20180901N/BDF/FLAT_0003_U.fits     0.2
/Users/oguzhan/tmp/20180901N/BDF/FLAT_0004_U.fits     0.2
/Users/oguzhan/tmp/20180901N/BDF/FLAT_0005_U.fits     0.2
/Users/oguzhan/tmp/20180901N/BDF/FLAT_0006_U.fits     0.2
/Users/oguzhan/tmp/20180901N/BDF/FLAT_0007_U.fits     0.2
/Users/oguzhan/tmp/20180901N/BDF/FLAT_0008_U.fits     0.2

```
