# -*- coding: utf-8 -*-
"""Atlantis astronomy framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

This is a set of function and classes to help with astronomical
calculations.

DISCLAIMER
    This software may be used, copied, or redistributed as long as it
    is not sold and this copyright notice is reproduced on each copy
    made. This routine is provided as is without any express or implied
    warranties whatsoever.

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    1 (2013-07-10 15:40 -0300)

REFERENCES
    Iris: A Python library for Meteorology and Climatology,
    http://scitools.org.uk/iris/

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

from numpy import (pi, meshgrid, deg2rad, arcsin, cos, sin, sign, asarray,
    recarray, arange, flatnonzero, rad2deg, round)

import sun

__all__ = ['Constants', 'Compass', 'sun', 'metergrid', 'gridcell_area']


class Constants:
    """Important geophysical constants in SI units.
    
    REFERENCES
        Moritz, H. Geodetic reference system 1980. Journal of Geodesy, 
        2000, 74, 128-162

    """
    # Earth's rotation rate, according to Moritz (2000)
    omega = 7292115e-11
    #omega = 2 * pi / (3600 * 23 + 56 * 60)
    days_in_year = 365.2421896698 # Wikipedia (?)
    hours_in_day = 2 * pi / omega / 3600
    minutes_in_day = hours_in_day * 60

    # Mean surface gravity [m / s**2], according to Moritz (2000)
    g = 9.797644656 # 9.81
    
    # Earth's mean radius [m], according to Moritz (2000). A spherical
    # Earth of radius 6371 km has the same volume as Earth. In this
    # case, one nautical mile is 1853.2488 m and one degree of latitude
    # is therefore 11194.93 m.
    a = 6371008.7714
    b = 111195.07973437 # Old 111177.5
    
    # Tidal frequencies [s**(-1)] in order of amplitude as in Stewart,
    # R. H. Introduction to physical oceanography; Texas A & M 
    # University, 2008, 345, available at http://oceanworld.tamu.edu/
    # resources/ocng_textbook/chapter17/chapter17_04.htm
    M2 = 2 * pi / (3600 * 12.4206);
    K1 = 2 * pi / (3600 * 23.9344);
    S2 = 2 * pi / (3600 * 12);
    O1 = 2 * pi / (3600 * 25.8194);


class Compass(object):
    """Compass class."""
    
    # 32 cardinal points of the compass
    # (source: http://en.wikipedia.org/wiki/Points_of_the_compass)
    _N_points = 32
    _n_points = 32
    _step = 1
    _delta = 11.25
    _cardinal_points = asarray(
        [
            ('North', 'N', 'Tramontana'),
            ('North by east', 'NbE', 'Qto Tramontana verso Greco'),
            ('North-northeast', 'NNE', 'Greco-Tramontana'),
            ('Northeast by north', 'NEbN', 'Qto Greco verso Tramontana'),
            ('Northeast', 'NE', 'Greco'),
            ('Northeast by east', 'NEbE', 'Qto Greco verso Levante'),
            ('East-northeast', 'ENE', 'Greco-Levante'),
            ('East by north', 'EbN', 'Qto Levante verso Greco'),
            ('East', 'E', 'Levante'),
            ('East by south', 'EbS', 'Qto Levante verso Scirocco'),
            ('East-southeast', 'ESE', 'Levante-Scirocco'),
            ('Southeast by east', 'SEbE', 'Qto Scirocco verso Levante'),
            ('Southeast', 'SE', 'Scirocco'),
            ('Southeast by south', 'SEbS', 'Qto Scirocco verso Ostro'),
            ('South-southeast', 'SSE', 'Ostro-Scirocco'),
            ('South by east', 'SbE', 'Qto Ostro verso Scirocco'),
            ('South', 'S', 'Ostro'),
            ('South by west', 'SbW', 'Qto Ostro verso Libeccio'),
            ('South-southwest', 'SSW', 'Ostro-Libeccio'),
            ('Southwest by south', 'SWbS', 'Qto Libeccio verso Ostro'),
            ('Southwest', 'SW', 'Libeccio'),
            ('Southwest by west', 'SWbW', 'Qto Libeccio verso Ponente'),
            ('West-southwest', 'WSW', 'Ponente-Libeccio'),
            ('West by south', 'WbS', 'Qto Ponente verso Libeccio'),
            ('West', 'W', 'Ponente'),
            ('West by north', 'WbN', 'Qto Ponente verso Maestro'),
            ('West-northwest', 'WNW', 'Maestro-Ponente'),
            ('Northwest by west', 'NWbW', 'Qto Maestro verso Ponente'),
            ('Northwest', 'NW', 'Maestro'),
            ('Northwest by north', 'NWbN', 'Qto Maestro verso Tramontana'),
            ('North-northwest', 'NNW', 'Maestro-Tramontana'),
            ('North by west', 'NbW', 'Qto Tramontana verso Maestro')
        ],
        dtype=[('compass_point', '|S28'), ('abbreviation', '|S28'),
            ('traditional_wind_point', '|S28')]
    )
    
    def __init__(self, point=0, N=32):
        """
        Initiates the compass class.

        PARAMETERS
            point (integer, float, string, optional) :
                Sets the current compass to the desired point.
                See `set_compass_point` for possible options.
            N (integer, optional) :
                Indicates the number of cardinal points in compass.
                Values can be either 4, 8, 16 or 32 (default).
        
        """
        if N in [4, 8, 16, 32]:
            self._n_points = int(N)
            self._step = 32 / self._n_points
            self._delta = 360. / self._n_points
            self._set_compass_point(point)
        else:
            raise ValueError('Invalid number of cardinal points.')


    def _set_compass_point(self, a, column=None):
        """
        Sets current compass point attributes.

        PARAMETER
            a (integer, float, string) :
                If input is an integer, assumes cardinal point
                number (1 to either 4, 8, 16 or 32). If input is a
                floating point number, then calculates the angle
                of the nearest cardinal point. If input is a string,
                searches in compass point abbreviation, compass point
                and traditional wind point names.
            column (string, array like, optional) :
                Column or list of columns to check for compass point.
        
        """
        if type(a) == int:
            self._n = a
        elif type(a) == float:
            # Makes sure that numers range from 0 to n-1.
            self._n = self.from_angle_to_point(a)
        elif type(a) in [str, unicode]:
            self._n = self._search_cardinal_point(a, column)
    
    
    def _search_cardinal_point(self, a, column=None):
        """Searches for cardinal point from string `a`."""
        if (column == None):
            _list = ['abbreviation', 'compass_point', 'traditional_wind_point']
        else:
            if type(column) in [list, tuple]:
                _list = column
            elif type(column) in [str, unicode]:
                _list = [column]
            else:
                raise ValueError('Invalid column `{}`.'.format(column))
        for _item in _list:
            try:
                return flatnonzero(
                    self._cardinal_points[_item][::self._step] == a
                )[0]
            except:
                pass
        #
        raise ValueError('Unable to locate compass point `{}`'.format(a))


    def _get_compass_point(self, column):
        return self._cardinal_points[column][::self._step][self._n]


    def _convert_angle(self, a, mode):
        """Converts angle `a` according to mode."""
        A, B = mode.split(':')
        if (A == 'degree' and a > 180):
            a -= 360
        if A == B:
            return a
        elif (A == 'degree') & (B == 'radian'):
            return deg2rad(a)
        else:
            raise ValueError('Invalid mode {}'.format(mode))


    def list_cardinal_points(self, n=None, mode='degree'):
        if n == None:
            n = arange(0, self._N_points, self._step)
            n = range(self._n_points)

        return [(i, self._cardinal_points[i*self._step][1],
            self._convert_angle(self._delta*i, mode='degree:{}'.format(mode)),
            self._cardinal_points[i*self._step][0],
            self._cardinal_points[i*self._step][2]) for i in n]


    def from_angle_to_point(self, a, mode='degree'):
        """
        Returns the compass point from input angle.

        PARAMETERS
            a (float, array like) :
                Single of list of angles
            mode (string, optional) :
                Determines the angle is in degree or radian.

        RETURNS
            index (integer, array like) :
                Single integer or list of integers giving the
                compass point.

        """
        if mode == 'radian':
            a = rad2deg(a)
        return (round(a / self._delta) % (self._n_points)).astype(int)
    

    @property
    def compass_point(self):
        """The compass point name."""
        return self._get_compass_point('compass_point')

    @compass_point.setter
    def compass_point(self, a):
        self._set_compass_point(str(a), 'compass_point')


    @property
    def abbreviation(self):
        """The compass point abbreviation."""
        return self._get_compass_point('abbreviation')

    @abbreviation.setter
    def abbreviation(self, a):
        self._set_compass_point(str(a), 'abbreviation')


    @property
    def traditional_wind_point(self):
        """The compass point traditional wind point."""
        return self._get_compass_point('traditional_wind_point')

    @traditional_wind_point.setter
    def traditional_wind_point(self, a):
        self._set_compass_point(str(a), 'traditional_wind_point')

    @property
    def angle(self):
        """The compass point angle in degrees."""
        return self._n * self._delta

    @angle.setter
    def angle(self, a):
        self._set_compass_point(float(a))


def metergrid(lon, lat, lon0=0, lat0=0, units='m'):
    """Converts zonal and meridional coordinates from degrees 
    latitude and longitude to another reference unit.
    
    PARAMETERS
        lon, lat (array like) :
            Longitude and latitude coordinates.
        lon0, lat0 (float, optional) :
            Reference longitude and latitude coordinates, default is 0.
        units (string, optional) :
            Units to which the coordinates will be converted. Accepted
            values are 'm' for meters (default), 'km' for kilometers,
            and 'nm' for International nautical miles.
    
    RETURNS
        x, y (array like) :
            New coordinates
    
    REFERENCES
        Haversine formula, available at
        http://en.wikipedia.org/wiki/Haversine_formula
    
    """
    # Makes sure longitude and latitude are numpy arrays
    lon = asarray(lon)
    lat = asarray(lat)
    # Checks if longitude and latitude are bi-dimensional
    if len(lon.shape) != len(lat.shape):
        raise ValueError('Number of dimensions of longitude and latitude '
                'arrays do not match.')
    if len(lon.shape) == 1:
        lon, lat = meshgrid(lon, lat)
    elif len(lon.shape) > 2:
        raise ValueError('Longitude and latitude arrays must be either one-'
                'dimensional or two-dimensional.')
    
    # Adapting the Haversine formula for regularly gridded arrays, the zonal
    # and meridional distances from the reference longitude \lambda_0 and
    # latitude \phi_0 are given by:
    #
    #   x = 2 r \arcsin{\left(\sqrt{cos^2(\phi) sin^2\left(\frac{\lambda -
    #       \lambda_0}{2}\right)}\right)}
    #   Y = 2 r \arcsin{\left(\sqrt{sin^2\left(\frac{\phi -
    #       \phi_0}{2}\right )}\right)}
    #
    # where r is the average radius of the spherical Earth, \lambda and \phi
    # are longitude and latitude in degrees.
    K = Constants()
    dlambda = deg2rad(lon - lon0)   # \Delta \lambda in radians.
    phi = deg2rad(lat)              # \phi in radians.
    dphi = deg2rad(lat - lat0)      # \Delta \phi in radians.
    x = K.a * cos(phi) * dlambda
    y = K.a * dphi
    if units == 'km':
        x *= 1e-3
        y *= 1e-3
    elif units == 'nm':
        x /= 1852
        y /= 1852
    return x, y


def gridcell_area(lon, lat):
    """Calculates the area of a grid cell whose center is given by
    longitude (lon) and latitude (lat).

    PARAMETERS
        lon, lat (array like):
            Longitude and latitude coordinates.

    RETURNS
        A (array like) :
            Grid cell area.
    
    """
    #
    # dA = r**2 * sin(\theta) * d{\theta} d{\phi}
    # \theta =
    # \phi = 
    #
    return False
