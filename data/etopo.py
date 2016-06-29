# -*- coding: utf-8 -*-
"""Atlantis data framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional data and
associated metadata. Much of this is based uppon Iris library

This module implements ETOPO1 Global Relief Model dataset reading
capabilities.

DISCLAIMER
    This software may be used, copied, or redistributed as long as it
    is not sold and this copyright notice is reproduced on each copy
    made. This routine is provided as is without any express or implied
    warranties whatsoever.

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    1 (2015-10-06 17:57 -0300)

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

from numpy import arange, argsort, flatnonzero, in1d, meshgrid
from scipy.io import netcdf_file as netcdf
from time import time

import atlantis.data
from klib.common import lon_n

DEBUG = True

class Grid(atlantis.data.Grid):
    """
    Common grid for General Bathymetric Chart of the Oceans (GEBCO).

    RESOLUTIONS
        30sec - 30 arc-seconds
        1min - 1 arc-minute
    
    """
    def __init__(self, path=None, resolution='1min', xlim=None, ylim=None):
        # Initializes the variables to default values. The indices 'n', 'k',
        # 'j' and 'i' refer to the temporal, height, meridional and zonal
        # coordinates respectively. If one of these indexes is set to 'None',
        # then it is assumed infinite size, which is relevant for the 'time'
        # coordinate.
        self.attributes = dict()
        self.dimensions = dict(n=0, k=0, j=0, i=0)
        self.coordinates = dict(n=None, k=None, j=None, i=None)
        self.variables = dict()
        self.params = dict()
        self.stencil_coeffs = dict()
        self.stencil_params = dict()
        
        # Sets global parameters for grid.
        if path == None:
            path = ('/media/academia/data/raw/etopo')
        self.params['path'] = path
        self.params['fname'] = 'ETOPO1_Bed_g_gmt4.grd'

        # Reads file header
        f = self.open_file()
        lon = f.variables['x'].data
        lat = f.variables['y'].data
        f.close()
        
        if (xlim != None) | (ylim != None):
            if xlim == None:
                xlim = (lon.min(), lon.max())
            if ylim == None:
                ylim = (lat.min(), lat.max())
            LON = lon_n(lon, xlim[0]) + 360
            i = argsort(LON)
            selx = i[flatnonzero((LON[i] >= xlim[0]) & (LON[i] <= xlim[1]))]
            sely = flatnonzero((lat >= ylim[0]) & (lat <= ylim[1]))
            selxx, selyy = meshgrid(selx, sely)
            lon = LON[selx]
            lat = lat[sely]
            self.params['xlim'] = xlim
            self.params['offset_i'] = selx[0]
            self.params['tesffo_i'] = selx[-1]
            self.params['selxx'] = selxx
            self.params['ylim'] = ylim
            self.params['offset_j'] = sely[0]
            self.params['tesffo_j'] = sely[-1]
            self.params['selyy'] = selyy

        # Sets the grid attributes, dimensions, coordinates and variables.
        self.name = 'ocean_bathymetry'
        self.title = f.title
        self.description = """
ETOPO1 is a 1 arc-minute global relief model of Earth's surface that
integrates land topography and ocean bathymetry. It was built from
numerous global and regional data sets, and is available in "Ice
Surface" (top of Antarctic and Greenland ice sheets) and "Bedrock" (base
of the ice sheets) versions.

Find out more about ETOPO at : http://www.ngdc.noaa.gov/mgg/global/global.html
"""
        self.dimensions = dict(n=1, k=1, j=lat.size, i=lon.size)
        self.coordinates = dict(n='time', k='height', j='latitude',
            i='longitude')
        self.variables = dict(
            time = atlantis.data.Variable(
                canonical_units='days since 0001-01-01 UTC',
                data=[0.],
            ),
            height = atlantis.data.get_standard_variable('height', data=[0.]),
            latitude = atlantis.data.get_standard_variable('latitude',
                data=lat),
            longitude = atlantis.data.get_standard_variable('longitude',
                data=lon),
            topo = atlantis.data.Variable(
                canonical_units='m',
                description='Topography of the Earth.'
            )
        )
        self.params['var_list'] = ['topo']
        
        return
    
    
    def create(self, *args, **kwargs):
        return


    def open_file(self):
        """Opens to netCDF file according to dataset parameters."""
        return netcdf('%s/%s' % (self.params['path'], self.params['fname']),
            'r')


    def read(self, t=None, z=None, y=None, x=None, N=None, K=None, J=None,
        I=None, var=None, nonan=True, result='full', profile=False,
        dummy=False):
        """Reads dataset.

        PARAMETERS
            t, z, y, x (array like, optional) :
                Sets the time, height, latitude and longitude for which
                the data will be read.
            N, K, J, I (array like, optional) :
                Sets the temporal, vertical, meridional and zonal
                indices for which the data will be read.
            var (string, optional) :
                Indicates which variable of the grid will be read. If
                the parameter is a list of variables, then the data will
                be returned as a list of arrays.
            nonan (boolean, optional) :
                If set to true (default) changes data values containing
                NaN to zero, preserving the mask.
            result (string, optional) :
                Determines wheter all time, height, latitude, longitude
                and data will be returned ('full', default), if
                temporal, vertical, meridional and zonal indices
                are returned instead ('indices'), or if only
                variable data is returned ('var only').
            profile (boolean, optional) :
                Sets whether the status is send to screen.
            dummy (boolean, optional) :
                If set to true, does not load data and returns the shape
                of the array that would have been returned.

        RETURNS
            t, z, y, x, dat (array like) :
                If 'result' is set to 'full', then all coordinates and
                data variables are returned.
            N, K, J, I, var (array like) :
                If 'result' is set to 'indices', then all indices and
                data variables are returned.
            dat (array like) :
                If 'result' is set to 'var only', then the data is
                returned.

        """
        global DEBUG
        t1 = time()
        
        # Checks input variables for consistency.
        if (t != None) & (N != None):
            raise ValueError('Both time and temporal index were provided.')
        if (z != None) & (K != None):
            raise ValueError('Both height and vertical index were provided.')
        if (y != None) & (J != None):
            raise ValueError(
                'Both latitude and meridional index were provided.')
        if (x != None) & (I != None):
            raise ValueError('Both latitude and zonal index were provided.')
        
        # Checks for variables indices. Intersects desired input values with
        # dataset dimesion data. In this dataset, since only surface data is
        # available, the height values are always zero.
        if t != None:
            N = flatnonzero(in1d(self.variables['time'].data, t))
        elif N == None:
            N = arange(self.dimensions['n'])
        if z != None:
            K = [0]
        elif K == None:
            K = [0]
        elif K != None:
            K = [0]
        if y != None:
            J = flatnonzero(in1d(self.variables['latitude'].data, y))
        elif J == None:
            J = arange(self.dimensions['j'])
        if x != None:
            I = flatnonzero(in1d(self.variables['longitude'].data, x))
        elif I == None:
            I = arange(self.dimensions['i'])

        # Sets the shape of the data array.
        shape = (len(N), 1, len(J), len(I))
        if dummy:
            return shape

        # Selects data according to indices.
        t = self.variables['time'].data[N]
        z = self.variables['height'].data
        y = self.variables['latitude'].data[J]
        x = self.variables['longitude'].data[I]
        xx, yy = meshgrid(x, y)
        II, JJ = meshgrid(I, J)

        try:
            sel_xx = self.params['selxx'][JJ, II]
            sel_yy = self.params['selyy'][JJ, II]
        except:
            sel_xx = II
            sel_yy = JJ
        
        # Read data
        f = self.open_file()
        var = f.variables['z'][sel_yy, sel_xx]
        f.close()
        
        if DEBUG:
            print 't: ', t
            print 'z: ', z
            print 'y:', y
            print 'x:', x
            print 'var: ', var
            print 'N: ', N
            print 'K: ', K
            print 'J: ', J
            print 'I:', I
            print 'shape: ', shape
        
        if result == 'full':
            return t, z, y, x, var
        elif result == 'indices':
            return N, K, J, I, var
        elif result == 'var only':
            return var
        else:
            raise Warning("Result parameter set imporperly to '%s', "
                "assuming 'var only'." % (result))
            return var
