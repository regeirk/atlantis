# -*- coding: utf-8 -*-
"""Rossby radius of deformation regridding.

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    1 (2013-07-05 15:51 -0300)

REFERENCES
    Chelton, D. B.; deSzoeke, R. A.; Schlax, M. G.; El Naggar, K. and
    Siwertz, N. Geographical variability of the first baroclinic Rossby
    radius of deformation. Journal of Physical Oceanography, American
    Meteorological Society, 1998, 28, 433-460.

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

import os
import numpy
import subprocess
from matplotlib import mlab
from scipy.io import netcdf_file as netcdf

import klib
from atlantis.data import dummy

# Parameters
icall = 'interpolate.sh'
nproc = 0

# Initializes dummy dataset
grid = dummy.Grid()
x = grid.variables['longitude']
y = grid.variables['latitude']

# Loads phase speed of first-mode baroclinic gravity waves and Rossby radius
# of deformation as of Chelton et al. (1998).
fname = 'rossrad.dat'
lat, lon, c, Ro = numpy.loadtxt(fname, unpack=True)
lon = klib.common.lon_n(lon, x.data.max())

# Saves data into temporary file
numpy.savetxt('dump_%d.xyz' % (nproc), numpy.array([lon, lat, Ro]).T)
ret = subprocess.call(['./%s' % (icall), '%s' % nproc], stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)

data = netcdf('%s_%d.grd' % ('output', nproc), 'r')
x = data.variables['x'].data
y = data.variables['y'].data
z = data.variables['z'].data
fname = 'rossrad_25.xy'
klib.file.save_map(x, y, z, '%s/%s.gz' % ('./', fname))

# Cleaning temporary files
os.remove('./dump_%d.xyz' % (nproc))
os.remove('./block_%d.xyz' % (nproc))
os.remove('./output_%d.grd' % (nproc))

# Plotting the results
crange = numpy.arange(0, 300, 25)
klib.gis.map(x, y, z, show=False, xlim=[x.min(), x.max()],
    units=r'\textnormal{km}', crange=crange, save='rossby_radius',
    extend='neither', orientation='worldmap', projection='moll')
