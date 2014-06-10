# -*- coding: utf-8 -*-
"""Rossby radius of deformation regridding.

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    2 (2013-11-04 18:32 -0300 DST)
    1 (2013-07-05 15:51 -0300)

REFERENCES
    Chelton, D. B.; deSzoeke, R. A.; Schlax, M. G.; El Naggar, K. and
    Siwertz, N. Geographical variability of the first baroclinic Rossby
    radius of deformation. Journal of Physical Oceanography, American
    Meteorological Society, 1998, 28, 433-460.

"""
from __future__ import division

__version__ = '$Revision: 2 $'
# $Source$

import klib

import os
import numpy
import subprocess
from matplotlib import mlab, pyplot
from scipy.io import netcdf_file as netcdf
from mpl_toolkits.basemap import cm

from atlantis.data import dummy
from atlantis.astronomy import metergrid

# Parameters
icall = 'interpolate.sh'
nproc = 0
etopo_file = '/home/sebastian/academia/data/etopo/etopo025.xy.gz'
mask_file = '/home/sebastian/academia/data/aviso/misc/mask025.xy.gz'
deg = numpy.array([0, 0.5, 1, 1.5, 2, 2.5])
deglbl = [r'0$^{\circ}$', r'0.5$^{\circ}$', r'1$^{\circ}$', r'1.5$^{\circ}$', 
    r'2$^{\circ}$', r'2.5$^{\circ}$']
#
pyplot.close('all')
pyplot.ion()

# Initializes dummy dataset
grid = dummy.Grid()
x = grid.variables['longitude'].data
y = grid.variables['latitude'].data

# Loads ETOPO 0.25 topography and mask
ex, ey, _, ez = klib.file.load_map(etopo_file, lon=x)
mx, my, _, mz = klib.file.load_map(mask_file, lon=ex, lat=ey)
xm, ym = metergrid([1.], ey, unit='km')

mz[numpy.isnan(mz)] = 0
mask = (mz == 0)

# Loads phase speed of first-mode baroclinic gravity waves and Rossby radius
# of deformation as of Chelton et al. (1998).
fname = 'rossrad.dat'
lat, lon, c, Ro = numpy.loadtxt(fname, unpack=True)
lon = klib.common.lon_n(lon, x.max())

# Saves data into temporary file and interpolates it using GMT (Smith & Wessel)
numpy.savetxt('dump_%d.xyz' % (nproc), numpy.array([lon, lat, Ro]).T)
ret = subprocess.call(['./%s' % (icall), '%s' % nproc], stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)

data = netcdf('%s_%d.grd' % ('output', nproc), 'r')
x = data.variables['x'].data
y = data.variables['y'].data
z = data.variables['z'].data
z = numpy.ma.masked_where(mask, z)
zx_deg = z / xm
zy_deg = z / (ym[1:]-ym[:-1]).mean()

fname = 'rossrad_25.xy'
klib.file.save_map(x, y, z, '%s/%s.gz' % ('./', fname))

# Cleaning temporary files
os.remove('./dump_%d.xyz' % (nproc))
os.remove('./block_%d.xyz' % (nproc))
os.remove('./output_%d.grd' % (nproc))

###############################################################################
# Plotting the results
###############################################################################
crange = numpy.arange(-8000, 9000, 1000)
klib.gis.map(ex, ey, ez, show=True, units=r'\textnormal{m}', crange=crange,
    save='etopo025', extend='both', orientation='worldmap', projection='moll',
    cmap=cm.GMT_relief, fillcontinents=False, etopo=False)

#crange = numpy.arange(0, 300, 25)
crange = numpy.array([0, 25, 50, 100, 200, 300])
klib.gis.map(x, y, z, show=True, xlim=[x.min(), x.max()],
    units=r'\textnormal{km}', crange=crange, save='rossby_radius', ftype='pdf',
    extend='neither', orientation='worldmap', projection='moll',
    ctype='contour', colors='k', cmap=None, cbar=False, fmt='%d')

crange = numpy.array([0, 0.3, 0.5, 1., 1.5, 2, 2.5, 3])
klib.gis.map(x, y, zx_deg, show=True, xlim=[x.min(), x.max()],
    units=r'^{\circ}', crange=crange, save='rossby_radius_deg', ftype='pdf',
    extend='neither', orientation='worldmap', projection='moll',
    ctype='contour', colors='k', cmap=None, cbar=False)

fig = klib.graphics.figure(fp=dict(figsize=[3.5, 4.5]), 
    ap=dict(left=0.23, bottom=0.17, right=0.93, top=0.89, wspace=0, hspace=0))
ax = klib.graphics.plot(z.mean(axis=1), y, xlabel='$L_{R,1}$',
    xunits=r'\textnormal{km}', yscale='deg', fig=fig, xlim=[0, 250], 
    ylim=[-90, 90], xtick='auto:5', ytick='auto:6', linewidth=2.0)
bx = ax.twiny()
bx.plot(zx_deg.mean(axis=1), y, 'k--', linewidth=1.5)
bx.set_xticks(deg)
bx.set_xticklabels(deglbl)
fig.savefig('rossby_radius_zonal.pdf', dpi=150)

# Thats all folks.
print 'Done.'
