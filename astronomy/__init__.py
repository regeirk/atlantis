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

from numpy import pi, meshgrid, deg2rad, arcsin, cos, sin, sign, asarray

import sun

__all__ = ['sun', 'constants', 'metergrid', 'gridcell_area']

class constants:
    """Important geophysical constants in SI units.
    
    REFERENCES
        Moritz, H. Geodetic reference system 1980. Journal of Geodesy, 
        2000, 74, 128-162

    """
    # Earth's rotation rate, according to Moritz (2000)
    omega = 7292115e-11
    #omega = 2 * pi / (3600 * 23 + 56 * 60)

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


def metergrid(lon, lat, lon0=0, lat0=0, unit='m'):
    """Converts zonal and meridional coordinates from degrees 
    latitude and longitude to another reference unit.
    
    PARAMETERS
        lon, lat (array like) :
            Longitude and latitude coordinates.
        lon0, lat0 (float, optional) :
            Reference longitude and latitude coordinates, default is 0.
        unit (string, optional) :
            Unit to which the coordinates will be converted. Accepted
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
    K = constants()
    dlambda = deg2rad(lon - lon0)   # \Delta \lambda in radians.
    phi = deg2rad(lat)              # \phi in radians.
    dphi = deg2rad(lat - lat0)      # \Delta \phi in radians.
    x = K.a * cos(phi) * dlambda
    y = K.a * dphi
    if unit == 'km':
        x *= 1e-3
        y *= 1e-3
    elif unit == 'nm':
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
