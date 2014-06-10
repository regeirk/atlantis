# -*- coding: utf-8 -*-
"""Atlantis data framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional data and
associated metadata. Much of this is based uppon Iris library

This module implements ocean color dataset reading capabilities.

DISCLAIMER
    This software may be used, copied, or redistributed as long as it
    is not sold and this copyright notice is reproduced on each copy
    made. This routine is provided as is without any express or implied
    warranties whatsoever.

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    1 (2013-07-28 00:51 -0300)

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

from bz2 import BZ2File
from matplotlib import dates
from numpy import (arange, argsort, array, asarray, flatnonzero, in1d, isnan, 
    ma, meshgrid, nan, ones)
from os import listdir, remove
from pyhdf.SD import SD, SDC
from string import atof
from sys import stdout
from time import time
from uuid import uuid1 as _uuid

import atlantis.data

from klib.common import basins, profiler, lon_n, reglist as _reglist
from atlantis.astronomy import metergrid

DEBUG = False

class Grid(atlantis.data.Grid):
    """Common grid for ocean optical properties fro different sensors.
    
        SENSORS
            SeaWiFS -- Sea-viewing Wide Field-of-view Sensor
            MODISA -- Moderate Resolution Imaging Spectroradiometer
    
        PROPERTIES (PRODUCTS)
            CHL / chla_a
        
        LEVEL / BINS
            l3m -- Level 3 mapped
    
    """
    def __init__(self, path=None, sensor='SeaWiFS', resolution='9km', 
        mask_file=None, xlim=None, ylim=None):
        # Initializes the variables to default values. The indices 'n', 'k', 'j'
        # and 'i' refer to the temporal, height, meridional and zonal coordinates
        # respectively. If one of these indexes is set to 'None', then it is
        # assumed infinite size, which is relevant for the 'time' coordinate.
        self.attributes = dict()
        self.dimensions = dict(n=0, k=0, j=0, i=0)
        self.coordinates = dict(n=None, k=None, j=None, i=None)
        self.variables = dict()
        self.params = dict()
        self.data = dict()
        self.stencil_coeffs = dict()
        self.stencil_params = dict()

        # Sets global parameters for grid.
        if path == None:
            path = '/home/sebastian/academia/data/oceancolor'
        self.params['path'] = '%s/%s' % (path, sensor)
        self.params['mask_file'] = mask_file
        self.params['uuid'] = str(_uuid())
        
        # Generates list of files, tries to match them to the pattern and to 
        # extract the time. To help understanding the naming convetion and 
        # pattern, see the following example:
        #   A20131612013168.L3m_8D_CHL_chlor_a_9km.bz2
        # resolution = '[0-9]+km'
        if sensor == 'SeaWiFS':
            sensor_prefix = 'S'
        elif sensor == 'MODISA':
            sensor_prefix = 'A'
        else:
            sensor = '.*'
        file_pattern = ('(%s)([0-9]{4})([0-9]{3})([0-9]{4})([0-9]{3}).(L3m)_'
            '(8D)_(CHL)_(chlor_a)_(%s).bz2') % (sensor_prefix, resolution)
        flist = listdir(self.params['path'])
        flist, match = _reglist(flist, file_pattern)
        self.params['file_list'] = flist

        # Reads first file in dataset to determine array geometry and 
        # dimenstions (lon, lat)
        HDF = self._open_HDF('%s/%s' % (self.params['path'], 
            self.params['file_list'][0]))
        HDF_att = HDF.attributes()
        lon = arange(HDF_att['Westernmost Longitude'], 
            HDF_att['Easternmost Longitude'], HDF_att['Longitude Step'])
        lat = arange(HDF_att['Northernmost Latitude'], 
            HDF_att['Southernmost Latitude'], -HDF_att['Latitude Step'])
        
        # If lon_0 is set, calculate how many indices have to be moved in 
        # order for latitude array to start at lon_0.
        if (xlim != None) | (ylim != None):
            if xlim == None:
                xlim = (lon.min(), lon.max())
            if ylim == None:
                ylim = (lat.min(), lat.max())
            #
            LON = lon_n(lon, xlim[0]) + 360
            i = argsort(LON)
            selx = i[flatnonzero((LON[i] >= xlim[0]) & (LON[i] <= xlim[1]))]
            sely = flatnonzero((lat >= ylim[0]) & (lat <= ylim[1]))
            ii, jj = meshgrid(selx, sely)
            lon = LON[selx]
            lat = lat[sely]
            self.params['xlim'] = xlim
            self.params['ylim'] = ylim
            self.params['lon_i'] = ii
            self.params['lat_j'] = jj
        
        # Creates a structured array for start year, start day, end year and 
        # end day. Aftwerwards, the dates are converted from julian day to 
        # matplotlib format, i.e. days since 0001-01-01 UTC.
        time_list = array([('%s-01-01' % (item[1]), atof(item[2]), 
            '%s-01-01' % (item[3]), atof(item[4])) for item in match], 
            dtype=[('start_year', 'a10'), ('start_day', 'f2'), 
            ('end_year', 'a10'), ('end_day', 'f2')])
        time_start = (dates.datestr2num(time_list['start_year']) + 
            time_list['start_day'] - 1)
        time_end = (dates.datestr2num(time_list['end_year']) + 
            time_list['end_day'] - 1)
        time_middle = 0.5 * (time_start + time_end)
        
        # Initializes the grid attributes, dimensions, coordinates and
        # variables.
        self.name = 'mass_concentration_of_chlorophyll_a_in_sea_water'
        self.description = ('Chlorophyll-a pigment concentration '
            'inferred from satellite visible light radiance measurements.')
        self.attributes['institution'] = HDF_att['Data Center']
        self.attributes['sensor name'] = HDF_att['Sensor Name']
        self.dimensions = dict(n=time_middle.size, k=0, j=lat.size, i=lon.size)
        self.coordinates = dict(n='time', k='height', j='latitude',
            i='longitude')
        self.variables = dict(
            time = atlantis.data.variable(),
            height = atlantis.data.get_standard_variable('height'),
            latitude = atlantis.data.get_standard_variable('latitude'),
            longitude = atlantis.data.get_standard_variable('longitude'),
            chla = atlantis.data.get_standard_variable(
                'mass_concentration_of_chlorophyll_a_in_sea_water'
            ),
            xm = atlantis.data.variable(),
            ym = atlantis.data.variable(),
        )
        self.variables['time'].data = time_middle
        self.variables['time'].canonical_units = 'days since 0001-01-01 UTC' 
        #
        self.variables['height'].data = 0.
        self.variables['latitude'].data = lat
        self.variables['longitude'].data = lon
        self.variables['chla'].canonical_units = 'mg m-3'
        #
        self.variables['xm'].canonical_units = 'km'
        self.variables['xm'].description = 'Zonal distance.'
        self.variables['ym'].canonical_units = 'km'
        self.variables['ym'].description = 'Meridional distance.'
        self.variables['xm'].data, self.variables['ym'].data = (
            metergrid(self.variables['longitude'].data, 
            self.variables['latitude'].data, unit='km')
        )
        return
    
    
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
            J = arange(self.dimensions['j'])
        if x != None:
            I = flatnonzero(in1d(self.variables['longitude'].data, y))
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
        var = ma.zeros(shape)
        # Walks through every time index and loads data range from maps.
        for n, T in enumerate(t):
            t2 = time()
            if profile:
                s = '\rLoading data... %s ' % (profiler(shape[0], n + 1, 0, 
                    t1, t2),)
                stdout.write(s)
                stdout.flush()
            # Uncompresses and reads HDF file
            HDF = self._open_HDF('%s/%s' % (self.params['path'], 
                self.params['file_list'][N[n]]))
            HDF_att = HDF.attributes()
            # Loads scientific dataset (SDS) and calculates the parameter 
            # value using scalling equation, slope and intercept.
            SDS_name = HDF.datasets().keys()[0]
            SDS = HDF.select(SDS_name)
            SDS_att = SDS.attributes()
            if (('lon_i' in self.params.keys()) &
                ('lat_j' in self.params.keys())):
                P = SDS[:, :][self.params['lat_j'], self.params['lon_i']][JJ,
                    II]
            else:
                P = SDS[:, :][JJ, II]
            P[P <= SDS_att['Fill']] = nan
            P = (SDS_att['Slope'] * P + SDS_att['Intercept'])
            P = ma.masked_where(isnan(P), P)
            if nonan:
                P.data[P.mask] = 0
            #
            var[n, 0, :, :] = P[None, None, :, :]
        
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
            return t, z, y, x, var
        elif result == 'indices':
            return N, K, J, I, var
        elif result == 'var only':
            return var
        else:
            raise Warning("Result parameter set imporperly to '%s', "
                "assuming 'var only'." % (result))
            return var
    
    
    def _open_HDF(self, file_path):
        """Opens compressed HDF file and returns pointer."""

        # STEP 1: Unzips image
        fileHDF = './%s_%s.hdf' % (self.params['uuid'], 'dump')
        fin = BZ2File(file_path, 'r')
        fout = open(fileHDF, 'wb')
        fout.write(fin.read())
        fin.close()
        fout.close()

        # STEP 2: Opens HDF file
        HDF = SD(fileHDF)
        
        # STEP 3: Cleanup
        remove(fileHDF)
        
        # LAST STEP: Retruns HDF pointer
        return HDF
