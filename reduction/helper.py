#!/usr/bin/env python

__all__ = ['FitsCollection']

from os import path
from glob import glob
from ccdproc import ImageFileCollection


# To make recursive file search in given directory.
class FitsCollection(ImageFileCollection):

    def __init__(self, **kwargs):
        if 'location' in kwargs:
            directory = path.join(kwargs['location'], '**', '*.*')

            filenames = glob(directory, recursive=True)
            kwargs['filenames'] = filenames

        super(FitsCollection, self).__init__(**kwargs)
