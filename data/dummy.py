# -*- coding: utf-8 -*-
"""Atlantis data framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional data and
associated metadata. Much of this is based uppon Iris library

This module implements a dummy dataset for test purposes only.

DISCLAIMER
    This software may be used, copied, or redistributed as long as it
    is not sold and this copyright notice is reproduced on each copy
    made. This routine is provided as is without any express or implied
    warranties whatsoever.

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    1 (2013-07-01 11:25 -0300)

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

from numpy import arange, meshgrid, nan, concatenate, sin, cos, deg2rad, sign, pi, loadtxt, round, argwhere, rad2deg
from numpy.ma import zeros, array, asarray, ones
from numpy.random import rand, randn
from matplotlib import dates
from os import path

import atlantis.data
from atlantis import astronomy
from atlantis.astronomy import sun
from klib import dynamics

DEBUG = False

class Grid(atlantis.data.Grid):
    def __init__(self, resolution=25, xlim=None, ylim=None):
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
        self._data = dict()
        self.stencil_coeffs = dict()
        self.stencil_params = dict()
        
        # Loads land / ocean mask
        _path = path.dirname(__file__)
        mask_file = '%s/mask%03d.xy.gz' % (_path, resolution)
        dat = loadtxt('%s' % (mask_file))
        lon = dat[0, 1:]
        lat = dat[1:, 0]
        mz = dat[1:, 1:]
        #lat = arange(-90., 90., dy)  + dy / 2.
        #lon = 20. + arange(0., 360., dx) - dx / 2.

        # If xlim and ylim are set, calculate how many indices have to be moved
        # in order for latitude array to start at xlim[0].
        lon, lat, xlim, ylim, ii, jj = self.getLongitudeLatitudeLimits(lon,
            lat, xlim, ylim)
        self.params['xlim'], self.params['ylim'] = xlim, ylim
        self.params['lon_i'], self.params['lat_j'] = ii, jj
        
        #
        self.name = 'sea_surface_height_above_sea_level'
        self.description = ('Dummy sea surface height anomaly data set, for '
            'test purposes only. The data set consists of a global 0.25 degree'
            ' grid with a seasonal cycle, westward propagating planetary wave '
            'field, meso-scale eddies and random noise. The temporal time-step'
            ' is one day.')
        self.dimensions = dict(n=None, k=0, j=lat.size, i=lon.size)
        self.coordinates = dict(n='time', k='height', j='latitude',
            i='longitude')
        self.variables = dict(
            time = atlantis.data.get_standard_variable('time', data=None),
            height = atlantis.data.get_standard_variable('height'),
            latitude = atlantis.data.get_standard_variable('latitude', 
                data=lat),
            longitude = atlantis.data.get_standard_variable('longitude', 
                data=lon),
            ssh = atlantis.data.get_standard_variable(
                'sea_surface_height_above_sea_level'),
            xm = atlantis.data.variable(),
            ym = atlantis.data.variable(),
            mask = atlantis.data.variable(
                description = ('Land and ocean mask defined as follows: land '
                    '(0), Atlantic Ocean (1), Pacific Ocean (2),  Indian Ocean'
                    ' (3), Caribbean Sea (4), Gulf of Mexico (5), Tasman Sea '
                    '(6), Bay of Bengal (7)'),
                data = mz[jj, ii]
            )
        )
        self.variables['height'].data = 0.
        #
        self.variables['xm'].canonical_units = 'm'
        self.variables['xm'].description = 'Zonal distance.'
        self.variables['ym'].canonical_units = 'm'
        self.variables['ym'].description = 'Meridional distance.'
        self.variables['xm'].data, self.variables['ym'].data = (
            astronomy.metergrid(self.variables['longitude'].data,
            self.variables['latitude'].data, units='km')
        )

        # Determines the function parameters for
        #
        # f(t) = a + b\,t + c\,\sin(\omega\,t + \phi) + 
        #   d\,\cos(\omega\,t + \phi)
        #
        i, j = self.dimensions['i'], self.dimensions['j']
        self.params['a'] = 100 * randn(j, i)
        self.params['b'] = 0.0005 + 0.001 * randn(j, i)
        self.params['c'] = 25 * rand(j, i)
        self.params['omega'] = 0.017202791695176148 + 0.001 * randn(j, i)
        self.params['phi'] = pi * rand(j, i)
    
    
    def read(self, t=None, z=None, y=None, x=None, N=None, K=None, J=None,
        I=None, var=None, result='full', profile=False, dummy=False,
        noise=True):
        """Reads dataset.

        PARAMETERS
            t, z, y, x (array like, optional) :
                Sets the time, height, latitude and longitude for which
                the data will be read.
            n, k, j, i (array like, optional) :
                Sets the temporal, vertical, meridional and zonal
                indices for which the data will be read.
            var (string, optional) :
                Indicates which variable of the grid will be read. If
                the parameter is a list of variables, then the data will
                be returned as a list of arrays.
            result (string, optional) :
                Determines wheter all time, height, latitude, longitude
                and data will be returned ('full', default), if
                temporal, vertical, meridional and zonal indices
                are returned instead ('indices'), or if only
                variable data is returned ('var only').
        
        RETURNS
            t, z, y, x, dat (array like) :
                If 'result' is set to 'full', then all coordinates anda
                data variables are returned.
            dat (array like) :
                If 'result' is set to 'var only', then the data is
                returned.

        """
        global DEBUG
        
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

        # Checks for variables indices.
        T0 = 693596.5
        if t != None:
            N = t - T0
        elif N == None:
            N = [0]
        if z != None:
            K = [0]
        elif K == None:
            K = [0]
        elif K != None:
            K = [0]
        if y != None:
            y = set(y)
            J = [j for j, item in enumerate(self.variables['latitude'].data)
                if item in y]
        elif J == None:
            J = arange(self.dimensions['j'])
        if x != None:
            x = set(x)
            I = [i for i, item in enumerate(self.variables['longitude'].data)
                if item in x]
        elif I == None:
            I = arange(self.dimensions['i'])

        # Sets the shape of the data array.
        shape = (len(N), 1, len(J), len(I))
        # Selects data according to indices.
        dt = 1.
        t = asarray(N) + T0
        z = self.variables['height'].data
        y = self.variables['latitude'].data[J]
        x = self.variables['longitude'].data[I]
        xx, yy = meshgrid(x, y)
        Im, Jm = meshgrid(I, J)
        xm = self.variables['xm'].data[Jm, Im]
        ym = self.variables['ym'].data[Jm, Im]
        var = zeros(shape)
        # Once the indices are all set, calculates the dummy data.
        a, b, c, omega, phi = (self.params['a'], self.params['b'],
            self.params['c'], self.params['omega'], self.params['phi'])
        var[:, 0, :, :] = self._function(t, a, b, c, omega, phi, t0=T0)
        if noise:
            var += self._noise(shape)

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


    def _function(self, t, a, b, c, omega, phi, t0=0):
        """
        Annual oscillation and trend function:
        
            f(t) = a + b\,t + c\,\sin(\omega\,t + \phi)
        
        PARAMETERS
            t -- Time.
            a -- Initial value.
            b -- Linear trend.
            c -- Oscillation aplitude.
            omega -- Angular frequency.
            phi -- Phase shift.
        
        """
        if len(t.shape) == 1:
            _t = t[:, None, None]
        else:
            raise Warning('Help!!! I don\'t know what to do next. :\'(')
        f = a + b * (_t - t0)  + c * sin(omega * _t + phi)
        #
        return f
    

    def _noise(self, shape, A=1.):
        """Random noise map.

        It creates a random noise map of given 'shape' from a normal
        distribution.
        
        """
        return A * randn(*shape)
