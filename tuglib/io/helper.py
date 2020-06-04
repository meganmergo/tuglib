#!/usr/bin/env python

__all__ = ['find_files']

import os
from glob import iglob


def find_files(location, extension=''):
    search_patterns = list()
    full_path = list()

    if location[-1] != os.sep:
        location += os.sep

    if isinstance(extension, (tuple, list)):
        for ext in extension:
            search_patterns.append(
                os.path.join(location, '**', '*.' + ext))
    else:
        if extension:
            search_patterns.append(
                os.path.join(location, '**', '*.' + extension))
        else:
            search_patterns.append(
                os.path.join(location, '**', '*'))

    for pattern in search_patterns:
        full_path += [filename for filename in iglob(pattern, recursive=True)
                      if os.path.isfile(filename)]

    return full_path
