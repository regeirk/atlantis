# -*- coding: utf-8 -*-
"""Atlantis units framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

This is a set of function and classes to help with unit conversions.

TODO
----
.. Read UDUNITS-2 documentation
    http://www.unidata.ucar.edu/software/udunits/udunits-2.2.19/doc/
    udunits/udunits2.html
.. Check https://code.google.com/p/cfunits-python/
.. Create an unit class for common usage;
.. Wrapper to convert more complicated units, i.e. `m s-1`, `W m-2`;
.. Unit aliases, i.e. `meter`, `metre`, `m`;
.. Unit multipliers, i.e. `milli`, `kilo`, `mega`, `hecto`;
.. UTF-8 comparison:
    http://stackoverflow.com/questions/3400171/python-utf-8-comparison

Disclaimer
----------
This software may be used, copied, or redistributed as long as it is not
sold and this copyright notice is reproduced on each copy made. This
routine is provided as is without any express or implied warranties
whatsoever.

Authors
-------
.. Sebastian Krieger (sebastian.krieger@usp.br)

Revision
--------
.. 1 (2014-09-26 16:21 -0300)

Reverences
----------
.. [1] Iris: A Python library for Meteorology and Climatology. Available
   at http://scitools.org.uk/iris/

"""
from __future__ import division, unicode_literals

__version__ = '$Revision: 1 $'
# $Source$

from numpy import (arcsin, arctan, asarray, cos, cosh, meshgrid, ndarray,
    polyval, sin, sinh, tan, tanh)
from numpy import ma

__all__ = ['fromAtoB']


def fromAtoB(data_in, A, B, vtype=None):
    """
    Converts from input data from unit A to unit B.

    Parameters
    ----------
    data_in : float, array like
        Input data to be converted.
    A : string
        Input data unit.
    B : string
        Desired output data unit.
    vtype : string, optional
        Determines the type of the variable. If type is `stdev`, no
        offset is applied, only scaling.

    Returns
    -------
    data_out : float, array like
        Converted data.

    """
    #
    #
    if A == B:
        return data_in
    elif (A == 'degC') & (B == 'K'):
        p = [1, 273.15]
    elif (A == 'K') & (B == 'degC'):
        p = [1, -273.15]
    elif (A == 'mbar') & (B == 'Pa'):
        p = [1e2, 0]
    elif (A == 'Pa') & (B == 'mbar'):
        p = [1e-2, 0]
    elif (A == 'hPa') & (B == 'Pa'):
        p = [1e2, 0]
    elif (A == 'Pa') & (B == 'hPa'):
        p = [1e-2, 0]
    elif (A in ['%', '1e-2']) & (B == '1'):
        p = [1e-2, 0]
    elif (A == '1') & (B in ['%', '1e-2']):
        p = [1e2, 0]
    elif (A == '1e-3') & (B == '1'):
        p = [1e-3, 0]
    elif (A == '1') & (B == '1e-3'):
        p = [1e3, 0]
    elif (A == 'ppm') & (B == '1'):
        p = [1e-6, 0]
    elif (A == '1') & (B == 'ppm'):
        p = [1e6, 0]
    elif (A == 'km h-1') & (B == 'm s-1'):
        p = [1./3.6, 0]
    elif (A == 'm s-1') & (B == 'km h-1'):
        p = [3.6, 0]
    elif (A == 'mm h-1') & (B == 'm s-1'):
        p = [1./3.6 * 1e-6, 0]
    elif (A == 'm s-1') & (B == 'mm h-1'):
        p = [3.6 * 1e6, 0]
    elif (A == 'm s-1') & (B == 'knot'):
        p = [1.9438444924574, 0]
    elif (A == 'knot') & (B == 'm s-1'):
        p = [0.51444444444, 0]
    elif (A in [u'µmol m-2 s-1', 'µmol m-2 s-1', '&#181;mol m-2 s-1']) & (B == 'mol m-2 s-1'):
        p = [1e-6, 0]
    elif (A == 'mol m-2 s-1') & (B == u'µmol m-2 s-1'):
        p = [1e6, 0]
    elif (A in [u'µg l-1', 'µg l-1', '&#181;g l-1']) & (B == 'kg m-3'):
        p = [1e-6, 0]
    elif (A == 'kg m-3') & (B == u'µg l-1'):
        p = [1e6, 0]
    elif (A == 'mg l-1') & (B == 'kg m-3'):
        p = [1e-3, 0]
    elif (A == 'kg m-3') & (B == 'mg l-1'):
        p = [1e3, 0]
    elif (A == 'mg m-3') & (B == 'kg m-3'):
        p = [1e-6, 0]
    elif (A == 'kg m-3') & (B == 'mg m-3'):
        p = [1e6, 0]
    elif (A == u'µmol m-3') & (B == 'mol m-3'):
        p = [1e-6, 0]
    elif (A == 'mol m-3') & (B == u'µmol m-3'):
        p = [1e6, 0]
    elif (A == u'µmol l-1') & (B == 'mol m-3'):
        p = [1e-9, 0]
    elif (A == 'mol m-3') & (B == u'µmol l-1'):
        p = [1e9, 0]
    # The following conversions are experimental!!!
    elif (A == 'rfu') & (B == '1'):
        p = [1, 0]
    elif (A == 'ppb') & (B == '1'):
        p = [1, 0]
    #elif (A == '') & (B == ''):
    else:
        raise ValueError('Unable to convert from `{}` to `{}`.'.format(A, B))
    #
    if vtype == 'stdev':
        p[1] = 0
    #
    if isinstance(data_in, ma.MaskedArray):
        return ma.array(polyval(p, data_in.data), mask=data_in.mask)
    else:
        return polyval(p, data_in)


def fromUTMtoLonLat(E, N, zone, hemisphere=1, datum='WGS84'):
    """
    Converts geographical units from Universal Transverse Mercator (UTM)
    conformal projection to longitude and latitude.

    Parameters
    ----------
    E, N : float, array like
        Easting and northing geographic Cartesial coordinates in meters.
    zone : integer
    hemisphere : char, integer, optional
        Either `N` or `+1` for northern hemisphere or `S` or `-1` for
        southern hemisphere.
    datum: string, optional

    Returns
    -------
    lon, lat : float, array like
        Longitude and latitude equivalent to UTM coordinates.
    k : float, array like
    gamma : float, array like

    References
    ----------
    .. [1] Universal Transverse Mercator coordinate system. Available at
       https://en.wikipedia.org/wiki/
       Universal_Transverse_Mercator_coordinate_system

    """
    # Checks for easting and northing parameter data type. If they are arrays
    # we have to create a meshgrid to perform all the calculations.
    if (isinstance(E, (ndarray, list, tuple)) |
        isinstance(N, (ndarray, list, tuple))):
        E, N = meshgrid(asarray(E), asarray(N))
    # Converts easting and northing parameters to kilometers.
    E, N = E * 1e-3, N * 1e-3

    # Constants and parameters
    a = 6378.137 # Equatorial radius of the earth in km.
    if hemisphere in ['N', 'n', 1]:
        N0 = 0
        hemisphere = 1.
    elif hemisphere in ['S', 's', -1]:
        N0 = 10000 # Wikipedia assumes km.
        hemisphere = -1.
    k0 = 0.9996
    E0 = 500  # Again, Wikipedia assumes km.
    f = 1./ 298357223536

    # Some calculated parameters
    n = f / (2 - f)
    A = a / (1 + n) * (1 + n**2/4 + n**4 / 64)
    beta = [1./2*n - 2./3*n**2 + 37./96/n**3, 1./48*n**2 + 1./15*n**3,
        17./480*n**3]
    delta = [2*n - 2./3*n**2 - 2*n**3, 7./3*n**2 - 8./5*n**3, 56./15*n**3]

    # Simple lambda functions
    beta_sincosh = lambda j, epsilon, eta: beta[j-1] * sin(2*j*epsilon) * cosh(2*j*eta)
    beta_cossinh = lambda j, epsilon, eta: beta[j-1] * cos(2*j*epsilon) * sinh(2*j*eta)
    beta_coscosh = lambda j, epsilon, eta: 2 * j * beta[j-1] * cos(2*j*epsilon) * cosh(2*j*eta)
    beta_sinsinh = lambda j, epsilon, eta: 2 * j * beta[j-1] * sin(2*j*epsilon) * sinh(2*j*eta)

    # Intermediate values
    epsilon = (N - N0) / (k0 - A)
    eta = (E - E0) / (k0 - A)
    epsilon_ = epsilon - (beta_sincosh(1, epsilon, eta) + beta_sincosh(2, epsilon, eta) + beta_sincosh(3, epsilon, eta))
    eta_ = eta - (beta_cossinh(1, epsilon, eta) + beta_cossinh(2, epsilon, eta) + beta_cossinh(3, epsilon, eta))
    sigma_ = 1 - (beta_coscosh(1, epsilon, eta) + beta_coscosh(2, epsilon, eta) + beta_coscosh(3, epsilon, eta))
    tau_ = (beta_sinsinh(1, epsilon, eta) + beta_sinsinh(2, epsilon, eta) + beta_sinsinh(3, epsilon, eta))
    chi = arcsin(sin(epsilon_) / cosh(eta_))

    # Finally
    phi = chi + (delta[0] * sin(2*1*chi) + delta[1] * sin(2*2*chi) +
        delta[2] * sin(2*3*chi)) # Latitude
    lambda0 = zone * 6. - 183. # Longitude of reference meridian
    lambda_ = lambda0 + arctan(sinh(eta_) / cos(epsilon_))
    k = k0 * A / a * ((1 + ((1 - n)/(1 + n) * tan(phi))**2) * (((cos(epsilon_))**2 + (sinh(eta_))**2) / (sigma_**2 + tau_**2)))**0.5
    gamma = hemisphere * arctan((tau_ + sigma_ * tan(epsilon_) * tanh(eta_)) / (sigma_ - tau_ * tan(epsilon_) * tanh(eta_)))
    #
    return lambda_, phi, k, gamma
