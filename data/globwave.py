# -*- coding: utf-8 -*-
"""Atlantis data framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional data and
associated metadata. Much of this is based uppon Iris library

This module implements GlobWave along-track (L2P) gridded altimetry
geophysical dataset.

DISCLAIMER
    This software may be used, copied, or redistributed as long as it
    is not sold and this copyright notice is reproduced on each copy
    made. This routine is provided as is without any express or implied
    warranties whatsoever.

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    1 (2014-05-28 19:54 -0300)

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

from matplotlib import dates
from numpy import (append, arange, array, asarray, ceil, dtype, flatnonzero, 
    float64, floor, hstack, isnan, ma, nan, ones, rec, sort, zeros)
from os import listdir, remove
from gzip import open as gzopen
from uuid import uuid1 as uuid
from scipy.io import netcdf_file as netcdf
from sys import stdout
from time import time

import json

import atlantis.data

from klib.common import profiler, reglist, lon360
from atlantis.astronomy import metergrid

DEBUG = False


class _NumpyAwareJSONEncoder(json.JSONEncoder):
    """."""
    def default(self, obj):
        if isinstance(obj, ndarray) and obj.ndim == 1:
            return [x for x in obj]
        else:
            return json.JSONEncoder.default(self, obj)


class Sequence(atlantis.data.Sequence):
    """Common sequencial data for Aviso along-track gridded delayed time
    products (sea level anomaly, absolute dynamic topography).
    
    """
    
    # Constant parameters
    _delays = dict(
        gdr = 'delayed-time',
        nrt = 'near-real-time'
    )
    _products = dict(
        altimeter = 'sea level anomaly',
        sar = 'synthetic aperture radar'
    )
    _missions = dict(
        altika = 'SARAL/AltiKa',
        cryosat2 = 'Cryosat-2',
        envisat = 'Envisat',
        ers1 = 'ERS-1',
        ers2 = 'ERS-2',
        geosat = 'Geosat',
        gfo = 'GFO',
        jason1 = 'Jason-1',
        jason2 = 'Jason-2',
        topex = 'TOPEX/Poseidon'
    )
    _labels = dict(
        sar = 'SAR',
        altimeter = 'ALT',
        altika = 'ALKA',
        cryosat2 = 'CRYO',
        envisat = 'ENVI',
        ers1 = 'ERS1',
        ers2 = 'ERS2',
        geosat = 'GEOS',
        gfo = 'GFO_',
        jason1 = 'JAS1',
        jason2 = 'JAS2',
        topex = 'TOPX'
    )

    # Variables
    params = None

    def __init__(self, level='l2p', product='altimeter', delay='gdr',
        missions=None, path=None, profile=True):
        """
        Initializes the dataset class for reading along-track gridded
        sequential data from the SSALTO/DUACS distributed by Aviso.

        PARAMETERS
            level (text, optional) :
                Process level, l2p (default) are preprocessed satellite
                wave data.
            product (text, optional) :
                Variable to be read (altimeter or sar)
            delay (text, optional) :
                Selects whether delayed time products (gdr, default)
                or near-real time (nrt) products are read.
            missions (text, array like, optional) :
                Determines the satellite missions to be selected (i. e.
                altika, cryosat2, envisat, ers1, ers2, geosat, gfo,
                jason1, jason2, topex) If set to 'none', all available
                missions are used.
            variable (text, optional) :
                Either 'vfec' for validated, filtered, sub-sampled and
                LWE-corrected; or 'vxxc' for validated, non-filtered,
                non-sub-sampled and LWE-corrected data.
            path (text, optional) :
                Path to the dataset files.
            profile (boolean, optional) :
                If true (default) shows status on screen.
        
        """
        t0 = time()
        # Checks all the input parameters for consistency
        if product not in self._products.keys():
            raise ValueError('Invalid product "%s".' % (product))
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
        
        # Initializes parameters and attributes in class variable
        self.attributes = dict()
        self.dimensions = dict(n=0, k=0, j=0, i=0)
        self.coordinates = dict(n=None, k=None, j=None, i=None)
        self.variables = dict()
        self.params = dict(
            level = level,
            product = product,
            delay = delay,
            missions = missions,
        )

        # Creates an universally unique identifiers (UUID) for this instance
        self.params['uuid'] = str(uuid())

        # Sets path and missing value parameters
        if path == None:
            path = '%s/%s/%s/%s' % ('/academia/data/raw/globwave', 
                level, product, delay)
        self.params['path'] = path
        self.params['missing_value'] = -9999.
        
        # Loads attributes from saved index
        #try:
        self.load_attributes()
        #except:
        #    self.attributes['time_mission'], self.attributes['time_dataset'] = self.make_index(profile=profile)
        #    self.save_attributes()
        
        # Updates dimensions, coordinates and creates time variable
        self.dimensions['n'] = len(self.attributes['time_dataset'])
        self.coordinates['n'] ='time'
        self.variables['time'] = atlantis.data.Variable(
            canonical_units = 'days since 0001-01-01 UTC',
            data = array(sorted(self.attributes['time_dataset'].keys())),
            height = atlantis.data.get_standard_variable('height', data=[0.]),
            latitude = atlantis.data.get_standard_variable('latitude'),
            longitude = atlantis.data.get_standard_variable('longitude'),
        )
        return None


    def load_attributes(self):
        """."""
        # Loads dataset attributes from JSON formatted file.
        fname = '%s/%s' % (self.params['path'], '.atlantis')
        with open(fname, 'r') as f:
            index_data = json.load(f)
        attribs = index_data['attributes']
        # Change the data type for some attributes' keys
        attribs['time_dataset'] = dict(
            (float(k), v) for k, v in attribs['time_dataset'].iteritems()
        )
        #
        self.attributes = attribs
        


    def save_attributes(self):
        """."""
        # Converts dataset information to JSON format and saves the data
        # to the dataset description file (.atlantis)
        fname = '%s/%s' % (self.params['path'], '.atlantis')
        dump = dict()
        dump['attributes'] = self.attributes
        url = '%s/%s' % (self.params['path'], '.atlantis')
        f = open(url, 'w')
        json.dump(dump, f, indent=2, cls=_NumpyAwareJSONEncoder)
        f.close()


    def make_index(self, profile=True):
        """."""
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
            tt1 = time()
            #
            mpath = '%s/%s' % (self.params['path'], mission)  # Mission path
            ylist = listdir(mpath)  # Year list in mission path
            Nyear = len(ylist)
            file_pattern = ('%s_%s_%s_%s_%s_(\d*)_(\d*)_(\d*)_(\d*)_(\d*)_'
                '(\d*).nc.gz') % ('GW', self.params['level'].upper(),
                self._labels[self.params['product']], self._labels[mission],
                self.params['delay'].upper())
            # Initializes time mission dictionary
            time_mission[mission] = dict(data=[], file=[])
            for j, yr in enumerate(ylist):
                tt2 = time()
                # Lists all the directories in year
                dlist = listdir('%s/%s' % (mpath, yr))
                for dset in dlist:
                    # Lists all the data files in mission in a given year and 
                    # matches it with the file pattern.
                    cur_path = '%s/%s/%s' % (mpath, yr, dset)
                    flist = listdir(cur_path)
                    flist.sort()
                    flist, match = reglist(flist, file_pattern)
                    # Convert data and product dates to matplotlib format, i.e. 
                    # days since 0001-01-01 UTC and appends to the global
                    # mission and dataset time dictionaries.
                    for k, item in enumerate(match):
                        datetime_start = dates.datestr2num(
                            '%4s-%2s-%2s %2s:%2s:%2s' % (item[0][0:4],
                            item[0][4:6], item[0][6:8], item[1][0:2],
                            item[1][2:4], item[1][4:6])
                        )
                        datetime_end = dates.datestr2num(
                            '%4s-%2s-%2s %2s:%2s:%2s' % (item[2][0:4],
                            item[2][4:6], item[2][6:8], item[3][0:2],
                            item[3][2:4], item[3][4:6])
                        )
                        time_data = (datetime_start + datetime_end) / 2.
                        cycle = int(item[4])
                        orbit = int(item[5])
                        time_mission[mission]['data'].append(time_data)
                        #
                        fname = '%s/%s/%s' % (yr, dset, flist[k])
                        descriptor = (mission, dset, fname, cycle, orbit)
                        if time_data not in time_dataset.keys():
                            time_dataset[time_data] = [descriptor]
                        else:
                            time_dataset[time_data].append(descriptor)
                        #
                        time_mission[mission]['file'].append(fname)
                #
                # Profiling
                if profile:
                    s = '\rBuilding preliminary time array for %s: %s ' % (
                        self._missions[mission], profiler(Nyear, j+1, t0, tt1,
                        tt2),
                    )
                    stdout.write(s)
                    stdout.flush()
            #
            time_mission[mission]['data'] = array(
                time_mission[mission]['data']
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

        return time_mission, time_dataset
    
    
    def read(self, var, x=None, y=None, radius=0., tlim=None, ylim=None,
        xlim=None, missions=None, sort=True, profile=True):
        """Reads dataset.
        
        PARAMETERS
            var (string) :
                Variable to be read from dataset. It also accepts
                special naming conventions in order to rename the
                original dataset variable and to load alternative
                variables in case of invalid data according to the
                syntax '[new_var_name]:var[|other_var]'.
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
        if var.find(':') >= 0:  # Checks spetial variable syntax
            var_name, var = var.split(':')
        else:
            var_name = var
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
            for (mission, dset, fname, cycle,
                orbit) in self.attributes['time_dataset'][tm]:
                # Skips mission not in missions list.
                if mission not in missions:
                    continue
                # Uncompresses gzipped file and opens NetCDF instance.
                data = self.read_file('%s/%s/%s' % (self.params['path'], 
                    mission, fname))
                # Reads variable from NetCDF file.
                raw_time = self.read_variable(data, 'time')
                raw_lat = self.read_variable(data, 'lat')
                raw_lon = self.read_variable(data, 'lon')
                raw_dat = self.read_variable(data, var)
                # Select relevant data range according to limit parameters
                sel_from_time = (
                    (raw_time >= min(tlim)) & (raw_time <= max(tlim))
                )
                if (ylim != None) | (xlim !=None):
                    sel_from_limits = ones(data.dimensions['time'], dtype=bool)
                else:
                    sel_from_limits = zeros(data.dimensions['time'],
                        dtype=bool)
                if ylim != None:
                    sel_from_limits = (sel_from_limits & 
                        ((raw_lat >= min(ylim)) & (raw_lat <= max(ylim))))
                if xlim != None:
                    sel_from_limits = (sel_from_limits & 
                        ((raw_lon >= min(xlim)) & (raw_lon <= max(xlim))))
                # Select relevant data according to points and search radius.
                sel_from_radius =  zeros(data.dimensions['time'], dtype=bool)
                for xx, yy in zip(x, y):
                    distance2 = ((raw_lat - yy) ** 2 + 
                        (raw_lon - lon360(xx)) ** 2)
                    sel_from_radius = sel_from_radius | (distance2 <= radius2)
                #
                sel_data = flatnonzero(sel_from_time & 
                    (sel_from_limits | sel_from_radius) & (~isnan(raw_dat)))
                _time = raw_time[sel_data]
                _lat = raw_lat[sel_data]
                _lon = raw_lon[sel_data]
                _dat = raw_dat[sel_data]
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
            ('longitude', float64), (var_name, float64), ('mission', '|S3')])
        
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


    def read_variable(self, data, var):
        """Reads variable data from NetCDF file pointer."""
        # Checks for special variable syntax. If the first variable has only
        # invalid values, uses the second, third, fourth, ...
        if var.find('|') >= 0:  # Checks spetial variable syntax
            var_list = var.split('|')
        else:
            var_list = [var]
        try:
            for item in var_list:
                dat = data.variables[item]
                if not (dat.data == dat._FillValue).all():
                    break
        except:
            pass
        # Determines if variable has a quality flag
        try:
            qf = data.variables[data.variables[var].quality_flag]
        except:
            qf = None
        # Loads bathymetry data
        H = data.variables['bathymetry']
        #
        dat0 = 0
        #
        try:
            scale_factor = dat.scale_factor
        except:
            scale_factor = 1.
        #
        if var == 'time':
            if dat.units == 'seconds since 1985-01-01 00:00:00.0':
                dat0 = dates.datestr2num('1985-01-01 00:00:00.0')
                scale_factor = 1. / 86400. # seconds per day!
            else:
                raise Warning("I don't know what to do!")
        #
        values = array(dat.data, dtype=float)
        # Masks invalid values and data over land
        try:
            values[(dat.data == dat._FillValue) | (H.data > 0)] = nan
        except:
            pass
        #
        return values * scale_factor + dat0
