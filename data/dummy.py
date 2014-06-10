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
from numpy.random import rand
from matplotlib import dates
from os import path

import atlantis.data
from atlantis import astronomy
from atlantis.astronomy import sun
from klib import dynamics

DEBUG = False

class Grid(atlantis.data.Grid):
    def __init__(self, dx=0.25, dy=0.25):
        # Initializes the variables to default values. The indices 'n', 'k', 'j'
        # and 'i' refer to the temporal, height, meridional and zonal coordinates
        # respectively. If one of these indexes is set to 'None', then it is
        # assumed infinite size, which is relevant for the 'time' coordinate.
        self.attributes = dict()
        self.dimensions = dict(n=0, k=0, j=0, i=0)
        self.coordinates = dict(n=None, k=None, j=None, i=None)
        self.variables = dict()
        self.params = dict()
        self._data = dict()
        self.stencil_coeffs = dict()
        self.stencil_params = dict()
        
        # Loads land / ocean mask
        mask_file = '/home/sebastian/academia/data/aviso/misc/mask025.xy.gz'
        dat = loadtxt('%s' % (mask_file))
        lon = dat[0, 1:]
        lat = dat[1:, 0]
        mz = dat[1:, 1:]
        #lat = arange(-90., 90., dy)  + dy / 2.
        #lon = 20. + arange(0., 360., dx) - dx / 2.
        
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
            time = None,
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
                data = mz
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
            self.variables['latitude'].data, unit='km')
        )
    
    
    def read(self, t=None, z=None, y=None, x=None, N=None, K=None, J=None,
        I=None, var=None, result='full', components=['seasonal', 'planetary',
        'eddy', 'noise']):
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
            components (list, optional) :
                A list containing which components will be included in
                the calculation. Options are the seasonal cycle
                ('seasonal'), westward propagating planetary waves
                ('planetary'), eddy fields ('eddy') and noise ('noise').

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
        for n, T in enumerate(t):
            d = zeros(shape[2:])
            if 'seasonal' in components:
                d += self._seasonal_cycle(T, xx, yy, A=4.)
            if 'planetary' in components:
                d += self._planetary_wave(T, xx, yy, A=[0.5, 1, 1., 0.5],
                    xm=xm, ym=ym)
            if 'eddies' in components:
                d += self._eddies(T, xx, yy, xm=xm, ym=ym)
            if 'noise' in components:
                d += self._noise(shape[2:])
            var[n, 0, :, :] = d[None, None, :, :]

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


    def _seasonal_cycle(self, t, x, y, A=1.):
        """Calculates a hypothetical seasonal cycle which is
        proportional to the astromical length of the day.

        """
        global DEBUG
        # Checks longitude and latitude dimensions.
        if x.shape != y.shape:
            raise ValueError('Longitude and latitude grid dimensions do not'
                ' match.')
        lod = length_of_day(t, x, y)
        
        if DEBUG:
            print 'lod: ', lod
        
        return A * (lod - 12.) / 12.


    def _planetary_wave(self, t, x, y, A=1., T=None, phase=None, xm=None,
        ym=None, mask=False):
        """Generates a hypothetical westward propagating planetary wave
        signal.

        PARAMETERS
            t (float) :
                Time in days.
            x, y (array like) :
                Longitude and latitude in degrees.
            A (float or array like, optional) :
                Amplitude of the Rossby wave signal. If given as a list,
                sets the amplitude for each wave period.
            T, phase (float or array like, optional) :
                Gives the wave period in days and wave phase in radians.
                If not set, assumes 3-, 6-, 12- and 24-month wave
                periods with zero phase.
            xm, ym (array like, optional) :
                Zonal and meridional grid in meters.
            mask (boolean, optional) :
                If true, masks equatorial region, where latitude is less
                than 5 degrees.
        
        """
        # Checks longitude and latitude dimensions.
        if x.shape != y.shape:
            raise ValueError('Longitude and latitude grid dimensions do not'
                ' match.')
        b, a = x.shape
        
        # Mask latitudes lower than 5 degrees.
        if mask:
            y = array(y, mask=(abs(y)<5.))
        
        # Calculate zonal and meridional grid in km, if necessary.
        if (xm == None) | (ym == None):
            xm, ym = astronomy.metergrid(x, y, unit='km')
        
        # Wave period in days and wave phase.
        if T == None:
            T = array([0.25, 0.5, 1., 2.]) * 365.25
        elif type(T) in [float, int]:
            T = array([T])
        elif type(T) in [list, tuple]:
            T = asarray(T)
        if phase == None:
            phase = arange(T.size)
        elif type(phase) in [float, int]:
            phase = array([phase])
        elif type(phase) in [list, tuple]:
            phase = asarray(phase)
        if T.shape != phase.shape:
            raise ValueError('Wave period and phase dimensions do not match.')
        if type(A) in [float, int]:
            A = array([A] * T.size)
        elif type(A) in [list, tuple]:
            A = asarray(A * int(T.size / len(A)))
        
        # Calculates theoretical phase speed [km day-1] for Rossby waves as
        # given in Oliveira & Polito (2013).
        phi = deg2rad(y)
        cp = -0.2 * abs(cos(phi)) / sin(phi)**2
        # cp = -0.2 * 1 / cos(phi)
        # Since cp = \frac{\omega}{k}, it is possible to construct wave like features
        # for the time periods T.
        Ro = zeros([b, a])
        for tt, pp, AA in zip(T, phase, A):
            omega = 2 * pi / tt
            k = omega / cp
            Ro += AA * cos(k * xm - omega * t)
        
        return Ro


    def _eddies(self, t, x, y, N=10, A=1., xm=None, ym=None):
        """Generates a hypothetical eddy field signal."""
        # Astronomical constants
        K = astronomy.constants()
        # Checks longitude and latitude dimensions.
        if x.shape != y.shape:
            raise ValueError('Longitude and latitude grid dimensions do not'
                ' match.')
        b, a = x.shape
        
        # Calculate zonal and meridional grid in km, if necessary.
        if (xm == None) | (ym == None):
            xm, ym = astronomy.metergrid(x, y, unit='km')
        
        # Loads first mode baroclinic Rossby radius of deformation file.
        _path = path.dirname(__file__)
        _fname = 'rossrad_25.xy.gz'
        dat = loadtxt('%s/%s' % (_path, _fname))
        Rx = dat[0, 1:]
        Ry = dat[1:, 0]
        RR = dat[1:, 1:]

        # Creates a couple of eddies with diameter given by the Rossby radius
        # of deformation and random amplitude.
        E = zeros([b, a])
        for n in range(N):
            j, i = round(rand(2) * [b, a])
            j, i = int(j), int(i)
            #
            dy, dx = abs(Ry - y[j, i]), abs(Rx - x[j, i])
            Ri = argwhere(dx == dx.min())[0][0]
            Rj = argwhere(dy == dy.min())[0][0]
            #
            dlambda = rad2deg(RR[Rj, Ri] / (K.a * cos(deg2rad(y[j, i]))))
            dphi = rad2deg(RR[Rj, Ri] / K.a)
            ii = argwhere((x[j, i] >= (x[j, i] - dlambda)) &
                (x[j, i] <= (x[j, i] + dlambda)))
            jj = argwhere((y[j, i] >= y[j, i] - dphi) &
                (y[j, i] <= y[j, i] + dphi))
            print i, j, x[j, i], y[j, i], RR[Rj, Ri], dlambda, dphi, ii, jj
            ii, jj = meshgrid(ii, jj)
            E[jj, ii] += A
            raise Warning('Continue from here!')
        
        return E

    def _noise(self, shape, A=1.):
        """Random noise map.

        It creates a random noise map of given 'shape' from a uniform
        distribution over [-0.5, 0.5).
        
        """
        return A * (rand(*shape) - 0.5)


def length_of_day(t, xx, yy):
    """Calculates the length of the day.

    The length of the day can be calculated as the time between sunrise
    and sunset which can be estimated from the sunrise equation. The
    length of the day is a function of latitude and time.

    PARAMETERS
        t (float) :
            Time in matplotlib date format (days since 0001).
        xx, yy (array like) :
            Longitude and latitude coordinate arrays (in degrees) as
            returned by numpy.meshgrid.

    RETURNS
        z (array like) :
            The duration of the day in hours for the given parameters.

    REFERENCES
        Sunrise Equation; Wikipedia; accessed on July 01, 2013; available
        at http://en.wikipedia.org/wiki/Sunrise_equation
    
    """
    S = sun.Sun()
    T =  dates.num2date(t)
    return S.dayLength(T.year, T.month, T.day, xx, yy)
