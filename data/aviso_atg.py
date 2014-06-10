# -*- coding: utf-8 -*-
"""Atlantis data framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional data and
associated metadata. Much of this is based uppon Iris library

This module implements Aviso along-track gridded sea level anomaly
dataset reading capabilities. It reads the gridded delayed time products
(sea level anomaly, geostrophic currents and error files).

DISCLAIMER
    This software may be used, copied, or redistributed as long as it
    is not sold and this copyright notice is reproduced on each copy
    made. This routine is provided as is without any express or implied
    warranties whatsoever.

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    1 (2014-05-15 12:48 -0300)

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

from matplotlib import dates
from numpy import (append, arange, array, asarray, ceil, dtype, flatnonzero, 
    float64, floor, hstack, rec, sort, zeros)
from os import listdir, remove
from gzip import open as gzopen
from uuid import uuid1 as uuid
from scipy.io import netcdf_file as netcdf
from sys import stdout
from time import time

import atlantis.data

from klib.common import profiler, reglist, lon360
from atlantis.astronomy import metergrid

DEBUG = False

class Sequence(atlantis.data.Sequence):
    """Common grid for Aviso along-track gridded delayed time products
    (sea level anomaly, absolute dynamic topography).
    
    """
    
    # Constant parameters
    _zones = {
        'global': 'global',
        'med': 'regional-mediterranean',
        'blacksea': 'regional-blacksea',
        'arctic': 'regional-arctic',
        'europe': 'regional-europe',
        'moz': 'regional-mozambique'
    }
    _delays = dict(
        dt = 'delayed-time',
        nrt = 'near-real-time'
    )
    _filterings = dict(
        vfec = 'filtered',
        vxxc = 'unfiltered'
    )
    _products = dict(
        sla = 'sea level anomaly',
        adt = 'absolute dynamic topography'
    )
    _missions = dict(
        e1 = 'ERS-1',
        e2 = 'ERS-2',
        tp = 'TOPEX/Poseidon',
        tpn = 'TOPEX/Poseidon new orbit',
        g2 = 'GFO',
        j1 = 'Jason-1',
        j1n = 'Jason-1 new orbit',
        j2 = 'Jason-2',
        en = 'Envisat',
        enn = 'Envisat new orbit',
        c2 = 'Cryosat-2',
        al = 'SARAL/AltiKa'
    )

    # Variables
    params = None

    def __init__(self, delay='dt', missions=None, zone='global',
        product='sla', variable='vxxc', path=None, profile=True):
        """
        Initializes the dataset class for reading along-track gridded
        sequential data from the SSALTO/DUACS distributed by Aviso.

        PARAMETERS
            delay (text, optional) :
                Selects whether delayed time products (dt, default)
                or near-real time products are read.
            missions (text, array like, optional) :
                Determines the satellite missions to be selected (i. e.
                e1, e2, tp, tpn, g2, j1, j1n, j2, en, enn, c2, al) If
                set to 'none', all available missions are used.
            zone (text, optional) :
                Geographic coverage of the selected products,
                    global -- Global geographic coverage;
                    med -- Mediterranean;
                    blacksea -- Black Sea;
                    moz -- Mozambique;
                    arctic -- Arctic;
                    europe -- Europe.
            product (text, optional) :
                Variable to be read (sla -- sea level anomaly or
                adt -- absolute dynamic topography)
            variable (text, optional) :
                Either 'vfec' for validated, filtered, sub-sampled and
                LWE-corrected; or 'vxxc' for validated, non-filtered,
                non-sub-sampled and LWE-corrected data.
            path (text, optional) :
                Path to the dataset files.
        
        """
        t0 = time()
        # Checks all the input parameters for consistency
        if delay not in self._delays.keys():
            raise ValueError('Invalid delay parameter "%s".' % (delay))
        if missions == None:
            missions = self._missions.keys()
        elif type(missions) == str:
            if missions in self._missions.keys():
                missions = [missions]
            else:
                raise ValueError('Invalid mission "%s".' % (missions))
        elif type(missions) == list:
            for item in missions:
                if item not in self._missions.keys():
                    raise ValueError('Invalid mission "%s".' % (item))
        else:
            raise ValueError('Invalid mission "%s".' % (missions))
        if zone not in self._zones.keys():
            raise ValueError('Invalid geographic zone "%s".' % (zone))
        if product not in self._products.keys():
            raise ValueError('Invalid product "%s".' % (product))
        if variable not in self._filterings.keys():
            raise ValueError('Invalid variable "%s".' % (variable))
        
        # Initializes parameters and attributes in class variable
        self.attributes = dict()
        self.dimensions = dict(n=0, k=0, j=0, i=0)
        self.coordinates = dict(n=None, k=None, j=None, i=None)
        self.variables = dict()
        self.params = dict(
            delay = delay,
            missions = missions,
            zone = zone,
            product = product,
            variable = variable
        )

        # Creates an universally unique identifiers (UUID) for this instance
        self.params['uuid'] = str(uuid())

        # Sets path and missing value parameters
        if path == None:
            path = '%s/%s/%s/%s/%s' % ('/academia/data/raw/aviso', 
                self._delays[delay], 'along-track', self._filterings[variable],
                product)
        self.params['path'] = path
        self.params['missing_value'] = -9999.
        
        # Determines the temporal range of the whole data set per mission
        t1 = time()
        if profile:
            s = '\rBuilding preliminary time array...'
            stdout.write(s)
            stdout.flush()
        
        time_mission = dict()
        time_dataset = dict()
        N = len(self.params['missions'])
        for i, mission in enumerate(self.params['missions']):
            t2 = time()
            #
            mpath = '%s/%s' % (path, mission)  # Mission path
            ylist = listdir(mpath)  # Year list in mission path
            file_pattern = '%s_%s_%s_%s_%s_(\d*)_(\d*).nc.gz' % (delay, zone, 
                mission, product, variable)
            time_mission[mission] = dict(data=[], product=[], file=[])
            for yr in ylist:
                # Lists all the data files in mission in a given year and 
                # matches it with the file pattern.
                flist = listdir('%s/%s' % (mpath, yr))
                flist.sort()
                flist, match = reglist(flist, file_pattern)
                # Convert data and product dates to matplotlib format, i.e. 
                # days since 0001-01-01 UTC and appends to the global mission
                # and dataset time dictionaries.
                for j, item in enumerate(match):
                    time_data = dates.datestr2num('%4s-%2s-%2s 12:00' % 
                        (item[0][:4], item[0][4:6], item[0][6:]))
                    time_mission[mission]['data'].append(time_data)
                    fname = '%s/%s' % (yr, flist[j])
                    descriptor = (mission, fname)
                    if time_data not in time_dataset.keys():
                        time_dataset[time_data] = [descriptor]
                    else:
                        time_dataset[time_data].append(descriptor)
                    #
                    time_product = dates.datestr2num('%4s-%2s-%2s 12:00' % 
                        (item[1][:4], item[1][4:6], item[1][6:]))
                    time_mission[mission]['product'].append(time_product)
                    #
                    time_mission[mission]['file'].append(fname)
            #
            time_mission[mission]['data'] = array(
                time_mission[mission]['data']
            )
            time_mission[mission]['product'] = array(
                time_mission[mission]['product']
            )
            time_mission[mission]['file'] = array(
                time_mission[mission]['file']
            )
            # Profiling
            if profile:
                s = '\rBuilding preliminary time array... %s ' % (profiler(N, 
                    i+1, t0, t1, t2),)
                stdout.write(s)
                stdout.flush()
        #
        if profile:
            stdout.write('\n')
            stdout.flush()
        #
        self.attributes['time_mission'] = time_mission
        self.attributes['time_dataset'] = time_dataset
        
        # Updates dimensions, coordinates and creates time variable
        self.dimensions['n'] = len(time_dataset)
        self.coordinates['n'] ='time'
        self.variables['time'] = atlantis.data.variable(
            canonical_units = 'days since 0001-01-01 UTC',
            data = array(sorted(time_dataset.keys())),
            height = atlantis.data.get_standard_variable('height', data=[0.]),
            latitude = atlantis.data.get_standard_variable('latitude'),
            longitude = atlantis.data.get_standard_variable('longitude'),
        )
        return None
    
    
    def read(self, x=None, y=None, radius=0., tlim=None, ylim=None, xlim=None,
            missions=None, sort=True, profile=True):
        """Reads dataset.
        
        PARAMETERS
            x, y (array like, optional) :
                List of zonal and meridional point coordinate of 
                interest.
            radius (float, optional) :
                Search radius in degrees.
            tlim, ylim, xlim  (array like, optional) :
                The temporal, meridional and zonal limits (minimum,
                maximum) for which data will be read.
            missions (array like, optional) :
                List of missions to read data from. If omitted, defaults
                available missions on dataset class intialization.
            sort (boolean optional) :
                If true, sorts the data record in order of ascendant 
                time, latitude and longitude.
            profile (boolean, optional) :
                Sets whether the status is send to screen.
        
        RETURNS
            dat (record array) :
                Record time-series of 'time', 'latitude', 'longitude', 
                selected variable and 'mission'.
        
        """
        t0 = time()
        # Checks input parameters.
        T = self.variables['time'].data
        if tlim == None:
            tlim = (T.min(), T.max())
        if (x != None) | (y != None):
            x, y = asarray(x), asarray(y)
            if x.size != y.size:
                raise ValueError('Zonal and meridional coordinate dimensions '
                    'do not match.')
            npoints = x.size
            radius2 = radius ** 2
        else:
            npoints = 0
            x = y = []
            #
            if ylim == None:
                ylim = (-90., 90.)
            if xlim == None:
                xlim = (0., 360.)
            else:
                # Make sure longitude limits are between 0 and 360.
                xlim = list(lon360(asarray(xlim)))
        if missions == None:
            missions = self.params['missions']
        
        # Aviso uses time in days since 1950-01-01 00:00:00 UTC, therefore
        # we have to calculate the initial time in matplotlib's format. We
        # also have to determine the proper variable using product name.
        T0 = dates.datestr2num('1950-01-01 00:00:00 UTC')
        var = self.params['product'].upper()
        
        # First we have to select which files will be loaded, which will 
        # depend on the temporal limits given in $t$.
        sel_time = flatnonzero((T >= floor(min(tlim))) & 
            (T <= ceil(max(tlim))))
        N = len(sel_time)
        
        # Second we will walk through each of the selected time in the dataset
        # and load the correspondant file for the available missions.
        t1 = time()
        if profile:
            s = '\rLoading data...'
            stdout.write(s)
            stdout.flush()
        # Reset important variables
        TIME, LAT, LON, VAR, MISSION = [array([])] * 5
        #
        for i, tm in enumerate(T[sel_time]):
            t2 = time()
            for (mission, fname) in self.attributes['time_dataset'][tm]:
                # Skips mission not in missions list.
                if mission not in missions:
                    continue
                # Uncompresses gzipped file and opens NetCDF instance.
                data = self.read_file('%s/%s/%s' % (self.params['path'], 
                    mission, fname))
                # Retrieve the scale factor for each variable
                scale_lat = data.variables['latitude'].scale_factor
                scale_lon = data.variables['latitude'].scale_factor
                scale_dat = data.variables[var].scale_factor
                # Get the raw time, latitude and longitude
                raw_time = data.variables['time'].data + T0
                raw_lat = data.variables['latitude'].data * scale_lat
                raw_lon = data.variables['longitude'].data * scale_lon
                # Select relevant data range according to limit parameters
                sel_from_time = (
                    (raw_time >= min(tlim)) & (raw_time <= max(tlim))
                )
                sel_from_limits = zeros(data.dimensions['time'], dtype=bool)
                if ylim != None:
                    sel_from_limits = (sel_from_limits | 
                        ((raw_lat >= min(ylim)) & (raw_lat <= max(ylim))))
                if xlim != None:
                    sel_from_limits = (sel_from_limits | 
                        ((raw_lon >= min(xlim)) & (raw_lon <= max(xlim))))
                # Select relevant data according to points and search radius.
                sel_from_radius =  zeros(data.dimensions['time'], dtype=bool)
                for xx, yy in zip(x, y):
                    distance2 = ((raw_lat - yy) ** 2 + 
                        (raw_lon - lon360(xx)) ** 2)
                    sel_from_radius = sel_from_radius | (distance2 <= radius2)
                #
                sel_data = flatnonzero(sel_from_time & 
                    (sel_from_limits | sel_from_radius))
                _time = raw_time[sel_data]
                _lat = raw_lat[sel_data]
                _lon = raw_lon[sel_data]
                _dat = data.variables[var].data[sel_data] * scale_dat
                #
                TIME = append(TIME, _time)
                LAT = append(LAT, _lat)
                LON = append(LON, _lon)
                VAR = append(VAR, _dat)
                MISSION = append(MISSION, [mission] * len(sel_data))
                #
                self.close_file(data)
            #
            # Profiling
            if profile:
                s = '\rLoading data... %s ' % (profiler(N, i+1, t0, t1, t2),)
                stdout.write(s)
                stdout.flush()
        #
        if profile:
            stdout.write('\n')
            stdout.flush()

        # Converts the data a structured array
        DAT = rec.fromarrays((TIME, LAT, LON, VAR, MISSION), 
            dtype=[('time', float64), ('latitude', float64), 
            ('longitude', float64), (self.params['product'], float64), 
            ('mission', '|S3')])
        #DAT = hstack((TIME[:, None], LAT[:, None], LON[:, None], 
        #    VAR[:, None], MISSION[:, None])).view(dtype=[('time', float64), 
        #    ('latitude', float64), ('longitude', float64), 
        #    (self.params['product'], float64), ('mission', '|S3')])
        
        # Some data sorting?
        if sort:
            DAT.sort(order=('time', 'latitude', 'longitude'), axis=0)
        
        return DAT


    def read_file(self, filename):
        """Reads zipped NetCDF file and returns its file pointer."""
        #
        # Uncompress NetCDF file.
        f = gzopen('%s' % (filename), 'rb')
        g = open('%s_%s.nc' % (self.params['uuid'], 'dump'), 'wb')
        g.write(f.read())
        f.close()
        g.close()
        #
        return netcdf('%s_%s.nc' % (self.params['uuid'], 'dump'), 'r')


    def close_file(self, data):
        """Closes and deletes temporary file."""
        #
        # Closes NetCDF file.
        try:
            data.close()
        except:
            pass
        
        # Removes the temporary dump file.
        try:
            remove('%s_%s.nc' % (self.params['uuid'], 'dump'))
        except:
            pass
