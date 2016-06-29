# -*- coding: utf-8 -*-
"""
Atlantis
========

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional data and
associated metadata. Much of this is based uppon Iris library

Authors
-------
.. Sebastian Krieger (sebastian.krieger@usp.br)

Revision history
----------------
.. 1 (2013-06-27 18:43 -0300)

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

import astronomy
import data
import dynamics
import units

__all__ = ['astronomy', 'data', 'dynamics', 'units']

# Documentation guide::
# https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt
