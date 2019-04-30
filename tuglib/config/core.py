#!/usr/bin/env python


class Config(object):

    def __init__(self):
        self.telescope = None
        self.instruments = None
        self.dedectors = None
        self.filters = None
        self.grisms = None
        self.apertures = None
        self.lamps = None

    def set_telescope(self, telescope):
        self.telescope = telescope

    def add_intrument(self, intrument):
        if self.instruments is None:
            self.instruments = list()

        self.instruments.append(intrument)

    def add_dedector(self, dedector):
        if self.dedectors is None:
            self.dedectors = list()

        self.dedectors.append(dedector)

    def add_filter(self, filter):
        if self.filters is None:
            self.filters = list()

        self.filters.append(filter)

    def add_grism(self, grism):
        if self.grisms is None:
            self.grisms = list()

        self.grisms.append(grism)

    def add_aperture(self, aperture):
        if self.apertures is None:
            self.apertures = list()

        self.apertures.append(aperture)

    def add_lamp(self, lamp):
        if self.lamps is None:
            self.lamps = list()

        self.lamps.append(lamp)


class Telescope(object):

    def __init__(self, name, f_number, focal_length):
        self.name = name
        self.f_number = f_number
        self.focal_length = focal_length


class Instrument(object):

    def __init__(self, name, imaging, spectra):
        self.name = name
        self.imaging = imaging
        self.spectra = spectra


class Detector(object):

    def __init__(self, name, shape, pixel_size):
        self.name = name
        self.shape = shape
        self.pixel_size = pixel_size


class Filter(object):

    def __init__(self, name, center_wavelength, band_pass):
        self.name = name
        self.center_wavelength = center_wavelength
        self.band_pass = band_pass


class Grism(object):

    def __init__(self):
        pass


class Aperture(object):

    def __init__(self, name, slit_type, shape):
        self.name = name
        self.slit_type = slit_type
        self.shape = shape


class Lamp(object):

    def __init__(self, name, center_wavelength, band_pass):
        self.name = name
        self.center_wavelength = center_wavelength
        self.band_pass = band_pass
