# -*- coding: utf-8 -*-
"""
Atlantis dynamics
=================

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

This is a set of function and classes to help with hydrodynamics.

TODO
----
. Mixing triangle
. Optimum Multiparameter (OMP) analysis
. Geostrophic velocity

Disclaimer
----------
This software may be used, copied, or redistributed as long as it is not
sold and this copyright notice is reproduced on each copy made. This
routine is provided as is without any express or implied warranties
whatsoever.

Authors
-------
.. Sebastian Krieger (sebastian.krieger@usp.br)

Revision history
----------------
.. 1 (2014-09-26 16:21 -0300)

References
----------
.. [1] Iris: A Python library for Meteorology and Climatology. Available
   at http://scitools.org.uk/iris/

"""
from __future__ import division, unicode_literals

__version__ = '$Revision: 1 $'
# $Source$

import gsw

from numpy import  asarray, c_, r_, ones, linalg, ma, interp

__all__ = ['water_mass_mixing', 'rho']


def brunt_vaisala(SP, t, p, lon=-45.401666, lat=-23.817233,  interpolate=True, 
    result='frequency'):
    """
    Returns Brunt-Väisäla (or buoyancy) frequency/period from practical
    salinity, in situ temperature and pressure (depth).
    
    Parameters
    ----------
    SP : array like
        Salinity (PSS-78) [1e-3]
    t : array like
        Temperature (ITS-90) [degC]
    p : array like
        Pressure [dbar]
    lon : array like, float, optional
        Longitude, decimal degrees east
    lat : array like, float, optional
        Latitude, decimal degrees north
    interpolate : boolean, optional
        If True, interpolates calculated frequency to original pressure
        (depth) levels. If False, returns results at mid pressure 
        between pressugre grid.
    result : {'freuency', 'period'}, optional
        Sets whether to return `frequency` [in 
    
    Returns
    -------
    N : array like
        Brunt-Väisälä frequency [s-1] or period [s].
    p : array like
        Pressure levels [dbar].
    
    
    """
    # Calculates Absolute Salinity from Practical Salinity
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    # Calcualtes Conservative Temperature from in situ temperature
    CT = gsw.CT_from_t(SA, t, p)
    #
    N2, p_mid = gsw.Nsquared(SA, CT, p, lat=lat)
    N2[N2 < 0] = 0
    #
    if interpolate:
        N2 = interp(p, p_mid, N2)
    #
    if result == 'frequency':
        return N2**0.5, p_mid
    elif result == 'period':
        return 1./N2**0.5, p_mid
    else:
        raise ValueError('Invalid result type {}.'.format(result))


def rho(SP, t, p, lon=-45.401666, lat=-23.817233):
    """
    Calculates in situ density from practical salinity, in situ 
    temperature and pressure (depth).
    
    Parameters
    ----------
    SP : array like
        Salinity (PSS-78) [1e-3]
    t : array like
        Temperature (ITS-90) [degC]
    p : array like
        Pressure [dbar]
    lon : array like, float, optional
        Longitude, decimal degrees east
    lat : array like, float, optional
        Latitude, decimal degrees north
    
    Returns
    -------
    rho : array like
        In situ density [kg m-3]
    
    """
    # Calculates Absolute Salinity from Practical Salinity
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    # Calcualtes Conservative Temperature from in situ temperature
    CT = gsw.CT_from_t(SA, t, p)
    # Calculates and returns density
    return gsw.rho(SA, CT, p)


def water_mass_mixing(T, S, indices):
    """
    Computes the water mass mixing percentage based on water mass core
    indices.

    The calculations are based on the mixing triagle (Mamayev 1975).

    Parameters
    ----------
    T : array like
        Temperature.
    S : array like
        Salinity.
    indices : array like
        A list/array with the core thermohaline indices for each water
        mass.

    Returns
    -------
    m1, m2, m3 : array like
        Relative composition for water masses 1, 2 and 3.

    Examples
    --------

    Reference
    ---------
    [1] Mamayev, O. I. (Ed.). Temperature -- Salinity Analysis fo World
        Ocean Waters Elsevier, 1975, 11.

    Notes
    -----
    This function is based upon code developed by Filipe Fernandes and
    available at https://ocefpaf.github.io/python4oceanographers/blog/
    2014/03/24/watermass/.

    """
    # Makes sure input parameters are numpy arrays
    T, S, indices = asarray(T), asarray(S), asarray(indices)

    # Creates linear system of three equations based on input parameters.
    a = r_[indices, ones((1, 3))]
    b = c_[T.ravel(), S.ravel(), ones(T.shape).ravel()].T
    m = linalg.solve(a, b)

    # The mixing indices 
    m1 = m[0].reshape(T.shape)
    m2 = m[1].reshape(T.shape)
    m3 = m[2].reshape(T.shape)

    # Mask values outside the mising triangle.
    m1 = ma.masked_outside(ma.masked_invalid(m1), 0, 1)
    m2 = ma.masked_outside(ma.masked_invalid(m2), 0, 1)
    m3 = ma.masked_outside(ma.masked_invalid(m3), 0, 1)

    return m1, m2, m3

