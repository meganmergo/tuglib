#!/usr/bin/env python

from .core import *

import warnings
from astropy.io.fits.verify import VerifyWarning
from astropy.wcs import FITSFixedWarning


# For fixing FITS header warning
warnings.simplefilter('ignore', category=VerifyWarning)
warnings.simplefilter('ignore', category=FITSFixedWarning)
