# -*- coding: utf-8 -*-
"""Atlantis data framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional data and
associated metadata. Much of this is based uppon Iris library

This module implements a framework for the distributed data archive
DADO running a MySQL database server.

DISCLAIMER
    This software may be used, copied, or redistributed as long as it
    is not sold and this copyright notice is reproduced on each copy
    made. This routine is provided as is without any express or implied
    warranties whatsoever.

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    1 (2013-07-28 12:19 -0300)

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

from time import time
t0 = time() # module start time for profiler

import MySQLdb

import atlantis.data

DEBUG = False

###############################################################################
# TODO
#
# createVariable, readVariable
# createLevel, readLevel
# createDataset, readDataset
# createProject, readProject
# createStation, readStation
# insertData, readData

###############################################################################
# PARAMETERS
#

###############################################################################
# CLASSES AND FUNCTIONS
#
class Grid(atlantis.data.Grid):
    def __init__(self, host='localhost', user='root', passwd='', db=''):
        self.db = MySQLdb.connect(host=host, user=user, passwd=passwd, db=db)
        self.cursor = self.db.cursor()
