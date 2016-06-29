# -*- coding: utf-8 -*-
"""Atlantis data framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional data and
associated metadata. Much of this is based uppon Iris library

This module implements common input and output (I/O) functions.

DISCLAIMER
    This software may be used, copied, or redistributed as long as it
    is not sold and this copyright notice is reproduced on each copy
    made. This routine is provided as is without any express or implied
    warranties whatsoever.

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    1 (2015-04-23 16:00 -0300)

"""
from __future__ import division
__version__ = '$Revision: 1 $'
# $Source$

from os import listdir
from os.path import isdir

__all__ = ['listdir_recursive']


def listdir_recursive(path, recursive=True):
    if isinstance(path, str):
        path = [path]

    # Ressets list of files in subdirectories
    file_list = []

    # Walks through every path and looks for files. If entry is a folder, adds
    # it to the list of paths, if apropriate.
    while(len(path) > 0):
        path_curr = u'{}'.format(path.pop(0).encode('utf-8'))
        files = listdir(path_curr)
        for f in files:
            url = u'{}/{}'.format(path_curr, f)
            if isdir(url):
                if recursive:
                    path.append(url)
            else:
                file_list.append(url)
    
    # Sorts the list of files.
    file_list.sort()
    # Done!
    return file_list
