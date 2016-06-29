# -*- coding: utf-8 -*-
"""Atlantis data framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional data and
associated metadata. Much of this is based uppon Iris library

This module implements Aviso FES tidal constituents dataset reading 
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
    1 (2015-12-09 09:00 +0200)

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

from matplotlib import dates
from numpy import (arange, argsort, array, asarray, flatnonzero, in1d, isnan, 
    ma, meshgrid, nan, ones, unique, hstack)
from os import listdir
from scipy.io import netcdf_file as netcdf
from string import atof
from sys import stdout
from time import time

import atlantis.data
from atlantis.astronomy import metergrid
from klib.common import num2ymd, reglist as reglist

DEBUG = False

class Grid(atlantis.data.Grid):
    """Common grid for Aviso FES tidal constituents."""
    def __init__(self, path=None, xlim=None, ylim=None):
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
        self.data = dict()
        self.stencil_coeffs = dict()
        self.stencil_params = dict()
        self.alias = dict()

        # Sets global parameters for grid.
        if path == None:
            path = ('/home/sebastian/academia/data/aviso/auxiliary/tide_model/'
                'fes2012_heights/fes2012/data')
        self.params['path'] = path
        self.params['var_list'] = []
        
        # Generates list of files, tries to match them to the pattern and to 
        # extract the time. To help understanding the naming convetion and 
        # pattern, see the following example:
        #   uwnd.2015.nc
        file_pattern = '(.*).([0-9]{4}).nc'
        flist = listdir(self.params['path'])
        flist, match = reglist(flist, file_pattern)
        self.params['file_list'] = flist
        
        # Gets list of variables from file match.
        _vars, _years = zip(*match)
        self.params['var_list'] = unique(_vars)
        self.params['year_list'] = unique(_years)
        
        # Loads data from first variable and loads longitude and latitude data.
        # We assume that all data is homogeneous throughout the dataset. Then 
        # walks through each year and loads time vector.
        _var = self.params['var_list'][0]
        for _i, _year in enumerate(self.params['year_list']):
            fname = '{}.{}.nc'.format(_var, _year)
            try:
                data.close()
            except:
                pass
            data = self._open_file(fname)
            #
            if _i == 0:
                lon = data.variables['lon'].data
                lat = data.variables['lat'].data
                time = data.variables['time'].data
            else:
                time = hstack([time, data.variables['time'].data])
        
        # Time in dataset is given in `hours since 1800-1-1 00:00:0.0` and we 
        # convert it to matplotlib's date format.
        if data.variables['time'].units == 'hours since 1800-1-1 00:00:0.0':
            self.params['t0'] = dates.date2num(dates.datetime.datetime(1800, 1, 1, 0, 0))
            time = self.params['t0'] + time / 24.
        
        # If lon_0 is set, calculate how many indices have to be moved in 
        # order for latitude array to start at lon_0.
        lon, lat, xlim, ylim, ii, jj = self.getLongitudeLatitudeLimits(lon,
            lat, xlim, ylim)
        
        self.params['xlim'], self.params['ylim'] = xlim, ylim
        self.params['lon_i'], self.params['lat_j'] = ii, jj
        
        # Initializes the grid attributes, dimensions, coordinates and
        # variables.
        self.name = 'ncep_reanalysis'
        self.description = ('NCEP Reanalysis project is analysis/forecast '
            'system to perform data assimilation using past data from 1979 '
            'owards.')
        self.attributes['institution'] = data.institution
        self.dimensions = dict(n=time.size, k=0, j=lat.size, i=lon.size)
        self.coordinates = dict(n='time', k='height', j='latitude',
            i='longitude')
        self.variables = dict(
            time = atlantis.data.Variable(),
            height = atlantis.data.get_standard_variable('height'),
            latitude = atlantis.data.get_standard_variable('latitude'),
            longitude = atlantis.data.get_standard_variable('longitude'),
            xm = atlantis.data.Variable(),
            ym = atlantis.data.Variable(),
        )
        #
        self.variables['time'].data = time
        self.variables['time'].canonical_units = 'days since 0001-01-01 UTC' 
        #
        self.variables['height'].data = 0.
        self.variables['latitude'].data = lat
        self.variables['longitude'].data = lon
        #
        self.variables['xm'].canonical_units = 'km'
        self.variables['xm'].description = 'Zonal distance.'
        self.variables['ym'].canonical_units = 'km'
        self.variables['ym'].description = 'Meridional distance.'
        self.variables['xm'].data, self.variables['ym'].data = (
            metergrid(self.variables['longitude'].data, 
            self.variables['latitude'].data, units='km')
        )
        #
        data.close()

        # Walks through each variable file for the first year, reads their 
        # attributes and adds to the dataset definition.
        _year = self.params['year_list'][0]
        for _var in self.params['var_list']:
            fname = '{}.{}.nc'.format(_var, _year)
            data = self._open_file(fname)
            for _key in data.variables.keys():
                if _key in ['time', 'level', 'lat', 'lon']:
                    continue
                try:
                    self.variables[_key] = atlantis.data.get_standard_variable(
                        data.variables[_key].standard_name, 
                        units=data.variables[_key].units
                    )
                except:
                    self.variables[_key] = atlantis.data.Variable(
                        units=data.variables[_key].units
                    )
                self.alias[_key] = _var
            data.close()
        #
        return
        
    def _open_file(self, fname):
        """
        Returns netCDF file according to dataset parameters and file 
        name.
        
        """
        return netcdf('{}/{}'.format(self.params['path'], fname), 'r')
    
    
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
            components (list, optional) :
                A list containing which components will be included in
                the calculation. Options are the seasonal cycle
                ('seasonal'), westward propagating planetary waves
                ('planetary'), eddy fields ('eddy') and noise ('noise').
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
            J = arange(self.dimensions['j'])[::-1]
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
        val = dict()
        
        # Gets list of years to analyse.
        YMD = num2ymd(t)
        years = unique(YMD[:, 0])
        
        # Gets list of variables to load.
        if var == None:
            var = list(set(self.variables.keys()) - 
                set(['time', 'height', 'latitude', 'longitude', 'ym', 'xm']))
            
        # Walks through every year loads data range from maps.
        I_var, I_year = len(var), len(years)
        _N = I_var * I_year
        for i_var, _var in enumerate(var):
            _val = ma.zeros(shape)
            #
            for i_year, _year in enumerate(years):
                t2 = time()
                if profile:
                    s = '\rLoading data... %s ' % (profiler(_N, 
                        i_var*I_var+i_year+1, 0, t1, t2),)
                    stdout.write(s)
                    stdout.flush()
                
                # Reads data file and extracts data values.
                fname = '{}.{:.0f}.nc'.format(self.alias[_var], _year)
                data = self._open_file(fname)
                
                # Reads time array from data and intersects it with time array
                # to be loaded. Use `in1d` to get the intersect indicies.
                _time = self.params['t0'] + data.variables['time'].data / 24.
                sel = flatnonzero(in1d(_time, t))
                les = flatnonzero(in1d(t, _time))
                                
                # Loads data and fits to desired grid.
                __var = data.variables[_var]
                
                if __var.data.ndim == 4:
                    __val = __var.data[sel, :, :, :]
                elif __var.data.ndim == 3:
                    __val = __var.data[sel, None, :, :]
                else:
                    raise ValueError('HELP! I crashed and I don\'t know what '
                        'to do next.')
                if ('lon_i' in self.params.keys()):
                    __val = __val[:, :, :, self.params['lon_i'][0, :][I]]
                else:
                    __val = __val[:, :, :, I]
                if ('lat_j' in self.params.keys()):
                    __val = __val[:, :, self.params['lat_j'][:, 0][J], :]
                else:
                    __val = __val[:, :, J, :]
                
                # Sets missing values to nan. But first we have to make sure 
                # that data in array are floating point numbers.
                __val = __val.astype(float)
                __val[(__val >= __var.missing_value)] = nan
                
                # And finally converts data according to scale factor and add 
                # offset. Please refer to the following website for further 
                # information 
                #   http://www.esrl.noaa.gov/psd/data/gridded/faq.html
                _val[les, :, :, :] = (__var.add_offset + 
                    __val * __var.scale_factor)

            # Masks invalid data.
            _val = ma.masked_invalid(_val)
            if nonan:
                _val.data[_val.mask] = 0
            #
            val[_var] = _val
        
        if len(var) == 1:
            val = val[var[0]]
        
        if profile:
            stdout.write('\r\n')
            stdout.flush()
        
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
            return t, z, y, x, val
        elif result == 'indices':
            return N, K, J, I, val
        elif result == 'var only':
            return val
        else:
            raise Warning("Result parameter set imporperly to '%s', "
                "assuming 'var only'." % (result))
            return var
