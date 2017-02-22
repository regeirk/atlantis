# -*- coding: utf-8 -*-
"""Atlantis data framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional gridded data,
multi-dimensional multi-point data, and associated metadata. Much of
this is based uppon Iris library.

The purpose of the CF conventions is to require conforming datasets to
contain sufficient metadata that they are self-describing in the sense
that each variable in the file has an associated description of what it
represents, including physical units if appropriate, and that each value
 can be located in space (relative to earth-based coordinates) and time.

Important terms and definitions from the OGC NetCDF Core adopted by the
present framework:
    attribute: Attributes  hold metadata. They contain information about
        properties of a variable or an entire dataset.
    coordinate: A coordinate variable is a variable that identifies a
        coordinate axis.
    data: The ISO definition of data is the reinterpretable
        representation of information in a formalized manner suitable
        for communication, interpretation, or processing
        [ISO/IEC 2382-1].
    data value: According to OGC Observations and Measurements, a data
        value is the result of a specialized event called an
        observation. An observation is an act of observing a property or
        phenomenon, with the goal of producing an estimate of the value
        of the property [OGC 07-022r1].
    dataset: According to ISO, a dataset is an identifiable collection
        of data [ISO 19101].
    dimension: Dimensions are used to specify variable shapes, common
        grids, and coordinate systems.
    global attribute (dataset attribute): Global attributes apply to a
        whole data set and may be used to record properties of all the
        data in a file, such as processing history or conventions used.
    grid: A grid is a network composed of two or more sets of curves in
        which the members of each set intersect the members of the other
        sets in an algorithmic way [ISO 19123].
    metadata: Metadata is data about data [ISO 19115].
    model: A model is an abstraction of some aspects of reality [ISO
        19109].
    record: A record is a finite, named collection of related items
        (objects or values).
    shape: The shape of a variable is specified with a list of zero or
        more dimensions.
    variable: A variable has a name, type, shape, attributes, and
        values.
    variable attribute: Variable attributes record the properties of one
        variable.

DISCLAIMER
    This software may be used, copied, or redistributed as long as it
    is not sold and this copyright notice is reproduced on each copy
    made. This routine is provided as is without any express or implied
    warranties whatsoever.

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    1 (2013-06-26 18:29 -0300)

REFERENCES
    Iris: A Python library for Meteorology and Climatology,
    http://scitools.org.uk/iris/

    Eaton, Brian, J. Gregory, B. Drach, K. Taylor, S. Hankin, J. Caron,
    R. Signell, P. Bentley, G. Rappa, H. Höck, A. Pamment, and M. Juckes;
    NetCDF Climate and Forecast (CF) Metadata Conventions, Version 1.6,
    5 December, 2011; http://cfconventions.org/Data/cf-conventions/
    cf-conventions-1.6/build/cf-conventions.html.

    Domenico, Ben and S. Nativi (Eds.); CF-netCDF3 Data Model Extension
    standard; Open Geospatial Consortium, Version 3.1, ref. OGC 11-165r2
    http://www.opengis.net/doc/is/netcdf-data-model-extension/1.0

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

from lxml import etree
from numpy import (argsort, array, zeros, arange, asarray, flatnonzero, isnan,
    nan, ndarray, savetxt, in1d, loadtxt, meshgrid, linalg, hstack, ones,
    bool_)
from numpy import ma
from scipy.misc import factorial
from os import path, listdir
from sys import stdout
from time import time
from warnings import warn

from klib.common import lon_n, profiler, num2ymd, reglist
from atlantis.astronomy import metergrid

import cf
import json

__all__ = ['cf', 'Grid', 'Variable']

DEBUG = False


# DATA FILE
# ~~~~~~~~~
#
# The general description of a file's contents should be contained in the
# following attributes: title, history, institution, source, comment and
# references.
#
#
# DIMENSIONS
# ~~~~~~~~~~~
#
#
# COORDINATES
# ~~~~~~~~~~~
#
# The use of Coordinate Variables is required for all dimensions that
# correspond to one dimensional space or time coordinates. All of a variable's
# dimensions that are latitude, longitude, vertical, or time dimensions must
# have corresponding coordinate variables, i.e., one-dimensional variables with
# the same name as the dimension.
#
# Coordinate Variable is defined as a numeric data type with values that are
# ordered monotonically. Missing values are not allowed in coordinate
# variables.
#
# 1. Latitude
#
#   float lat(lat) ;
#       lat:long_name = "latitude" ;
#       lat:units = "degrees_north" ;
#       lat:standard_name = "latitude" ;
#
#
# 2. Longitude
#
#   float lon(lon) ;
#       lon:long_name = "longitude" ;
#       lon:units = "degrees_east" ;
#       lon:standard_name = "longitude" ;
#
#
# 3. Vertical (height, depth, dimensional)
#
#   float lev(lev) ;
#       lev:long_name = "sigma at layer midpoints" ;
#       lev:positive = "down" ;
#       lev:standard_name = "atmosphere_sigma_coordinate" ;
#       lev:formula_terms = "sigma: lev ps: PS ptop: PTOP" ;
#
#   float pres(pres) ;
#       pres:long_name = "pres
#       pres:units = "hPa" ;
#
#
# 4. Time
#
#   double time(time) ;
#       time:long_name = "time" ;
#       time:units = "days since 1990-1-1 0:0:0" ;
#       time:calendar = "none" ;
#
#
# VARIABLES
# ~~~~~~~~~
#
# Each variable has an associated description which is provided by the
# attributes units, long_name, and standard_name.
#
#   units: The units attribute is required for all variables that represent
#       dimensional quantities (except for boundary variables defined. The
#       values of the units attributes are character strings that are
#       recognized by UNIDATA's Udunits package [UDUNITS].
#   long_name, standard_name: These attributes are used to describe the content
#       of each variable. Use of at least one of them is strongly recommended.
#       The use of standard names will facilitate the exchange of climate and
#       forecast data by providing unambiguous identification of variables most
#       commonly analyzed.
#
#
# BRAINSTORMING
# ~~~~~~~~~~~~~
# (...) mode parameter which must be a string with the name of the boundary
# condition. Following boundary conditions are currently supported: 'nearest',
# 'wrap', 'reflect', 'constant'
# (http://docs.scipy.org/doc/scipy/reference/tutorial/ndimage.html)
#


class Variable(object):
    # Some initializations
    _attr_list = ['long_name', 'standard_name', 'units', 'canonical_units',
        'description', 'grib', 'amip', 'id', 'missing_value', 'scale_factor',
        'add_offset', 'grids', 'data', 'standard_deviation', 'flag_values',
        'flag_meanings', 'string_format', 'symbol', 'valid_range']


    def __init__(self, **kwargs):
        """
        Attributes
        ----------
        long_name : string
            A long descriptive name which may, for example, be used for
            labeling plots.
        standard_name : string
            Provides unambiguous identification of variables most
            commonly analyzed.
        units : string
            A character string that specifies the units used for the
            variable's data.
        canonical_units : string
            Representative units of the physical quantity. These are the
            units in which data is stored.
        description : string
            The description is meant to clarify the qualifiers of the
            fundamental quantities.
        grib :
            GRIB parameter codes from ECMWF and NCEP.
        amip :
            AMIP identifier.
        id : int, optional
            Identification number of the variable. This is usefull when
            using databases.
        missing_value : float, complex
            Can be a scalar or vector (complex) containing values
            indicating missing data.
        scale_factor: float, optional
            Factor by which data should be multiplied after being read.
        add_offset: float, optional
            Offset to be added to data after being read. If both
            scale_factor and add_offset are given, data is first scaled.
        grids :
        data : array like
        flag_values : array like
            Contains a list of the possible flag values.
        flag_meanings : string
            A blank separated list of descriptive words or phrases, one
            for each flag value.

        """
        #
        for key, item in kwargs.items():
            try:
                setattr(self, key, item)
            except:
                print 'Warning: Invalid attribute {0}'.format(key)
                pass
        # Initializes other attributes if not given in kwargs.
        for key in set(self._attr_list) - set(kwargs.keys()):
            setattr(self, key, None)
        #
        return


    def __getattr__(self, name):
        # Checks whether attribute name is valid:
        if name not in self._attr_list:
            raise ValueError('Invalid attribute `{}`.'.format(name))
        # Returns attribute value or None if value was not assigned.
        try:
            return self.__dict__[name]
        except:
            return None


    def __setattr__(self, name, value):
        # Checks whether attribute name is valid:
        if name not in self._attr_list:
            raise ValueError('Invalid attribute `{}`.'.format(name))
        # Assigns value to class' dict of attributes.
        self.__dict__[name] = value


    @property
    def long_name(self):
        """The name of the variable."""
        return self.__getattr__('long_name')

    @long_name.setter
    def long_name(self, value):
        self.__setattr__('long_name', value)

    @property
    def standard_name(self):
        """The standard name of the variable."""
        return self.__getattr__('standard_name')

    @standard_name.setter
    def standard_name(self, value):
        self.__setattr__('standard_name', value)


    @property
    def units(self):
        """The units of the variable."""
        return self.__getattr__('units')

    @units.setter
    def units(self, value):
        self.__setattr__('units', value)


    @property
    def canonical_units(self):
        """The canonical_units of the variable."""
        return self.__getattr__('canonical_units')

    @canonical_units.setter
    def canonical_units(self, value):
        self.__setattr__('canonical_units', value)

    @property
    def description(self):
        return self.__getattr__('description')

    @description.setter
    def description(self, value):
        self.__setattr__('description', value)


    @property
    def grib(self):
        return self.__getattr__('grib')

    @grib.setter
    def grib(self, value):
        self.__setattr__('grib', value)


    @property
    def amip(self):
        return self.__getattr__('amip')

    @amip.setter
    def amip(self, value):
        self.__setattr__('amip', value)


    @property
    def missing_value(self):
        return self.__getattr__('missing_value')

    @missing_value.setter
    def missing_value(self, value):
        self.__setattr__('missing_value', value)


    @property
    def valid_range(self):
        return self.__getattr__('valid_range')

    @valid_range.setter
    def valid_range(self, value):
        self.__setattr__('valid_range', value)


    @property
    def scale_factor(self):
        """
        The scale factor by which data should be multiplied after being
        read.

        """
        return self.__getattr__('scale_factor')

    @scale_factor.setter
    def scale_factor(self, value):
        self.__setattr__('scale_factor', value)


    @property
    def add_offset(self):
        """
        The offset to be added to the data after being read. If both
        scale_factor and add_offset attributes are present, data is
        first scaled and then offset.

        """
        return self.__getattr__('add_offset')

    @add_offset.setter
    def add_offset(self, value):
        self.__setattr__('add_offset', value)


    @property
    def data(self):
        return self.__getattr__('data')

    @data.setter
    def data(self, value):
        self.__setattr__('data', value)


    @property
    def grids(self):
        """The grid description of the variable."""
        return self.__getattr__('grids')

    @grids.setter
    def grids(self, value):
        self.__setattr__('grids', value)


    # Extends default atlantis.data.variable class to include
    # database indexes and variable symbol.
    @property
    def id(self):
        """The ID number of the variable."""
        return self.__getattr__('id')

    @id.setter
    def id(self, value):
        self.__setattr__('id', value)


    @property
    def symbol(self):
        """The symbol of the variable."""
        return self.__getattr__('symbol')

    @symbol.setter
    def symbol(self, value):
        self.__setattr__('symbol', value)


    @property
    def string_format(self):
        """Defines format for data output."""
        return self.__getattr__('string_format')

    @string_format.setter
    def string_format(self, value):
        self.__setattr__('string_format', value)


    def json_dict(self):
        """
        Returns the variable representation as a dictionary for JSON...

        """
        dump = dict()
        for key, value in self.__dict__.items():
            if type(value) == ndarray:
                dump[key] = list(value)
            else:
                dump[key] = value
        return dump


class Sequence(object):
    """Defines a sequential dataset."""
    def __init__(self, path=None, pattern=None):
        self.params = dict(path=path, pattern=pattern)


    def read(self):
        """Reads sequence."""
        return False


    def write(self):
        """Writes sequence."""
        return False


class Grid(object):
    """Defines a spatio temporal dataset grid."""

    def _message(self, s):
        """Writes message to standard output."""
        stdout.write(s)
        stdout.flush()


    def __init__(self, path=None, pattern=None, xlim=None, ylim=None):
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
        self.default_vars = ['time', 'height', 'latitude', 'longitude', 'xm',
            'ym']

        if path == None:
            return
        else:
            self.params['path'] = path

        # Tries to read the dataset description file (.atlantis) to initialize
        # the variables.
        try:
            url = '%s/%s' % (self.params['path'], '.atlantis')
            f = open(url, 'r')
            d = json.load(f)
            f.close()
            loaded = True
        except:
            loaded = False
            pass

        if loaded:
            try:
                self.name = d['name']
                self.description = d['description']
                self.attributes = d['attributes']
                self.dimensions = d['dimensions']
                self.coordinates = d['coordinates']
                for var, params in d['variables'].items():
                    self.variables[var] = Variable()
                    self.variables[var].__dict__ = params
            except:
                raise IOError('Unable to read dataset description file.')

            for key in self.variables.keys():
                try:
                    if type(self.variables[key].data) == list:
                        self.variables[key].data = asarray(
                            self.variables[key].data
                        )
                except:
                    pass

            # If lon_0 is set, calculate how many indices have to be moved in
            # order for latitude array to start at lon_0.
            lon = self.variables['longitude'].data
            lat = self.variables['latitude'].data
            lon, lat, xlim, ylim, ii, jj = self.getLongitudeLatitudeLimits(lon,
                lat, xlim, ylim)
            self.params['xlim'], self.params['ylim'] = xlim, ylim
            self.params['lon_i'], self.params['lat_j'] = ii, jj
            self.variables['longitude'].data = lon
            self.variables['latitude'].data = lat
            self.dimensions['j'] = lat.size
            self.dimensions['i'] = lon.size

            # Create longitude and latitude grid in km.
            self.createMeterGrid()

        # List files in path matching pattern. The default pattern is
        # '(var)_([0-9]*).xy.gz', where 'var' is the short name of the
        # variable, i.e. 'ssh', 'sst', 'chla', etc.
        if loaded:
            var_list = list(set(self.variables.keys()) -
                set(self.default_vars))
            self.params['var_list'] = var_list
            if pattern == None:
                pattern = '(%s)_([0-9]*).xy.gz' % ('|'.join(var_list))
            self.params['pattern'] = pattern
            #
            n_var = len(var_list)
            for var in var_list:
                if n_var > 1:
                    fpath = '{0}/{1}'.format(self.params['path'], var)
                    rel_path = '{0}/'.format(var)
                else:
                    fpath = self.params['path']
                    rel_path = ''
                flist = listdir(fpath)
                flist.sort()
                flist, match = reglist(flist, self.params['pattern'])
                self.params['file_list_{}'.format(var)] = [
                    '{0}{1}'.format(rel_path, item) for item in flist
                ]
                if len(flist) != self.dimensions['n']:
                    warn(('List of files in \'{}\' does not match temporal '
                        'dimension length.').format(var))

        return


    def create(self, var=[]):
        """Creates the Atlantis dataset."""
        # Converts dataset information to JSON format and saves the data
        # to the dataset description file (.atlantis)
        fname = '%s/%s' % (self.params['path'], '.atlantis')
        dump = dict()
        dump['name'] = self.name
        dump['description'] = self.description
        dump['attributes'] = self.attributes
        dump['dimensions'] = self.dimensions
        dump['coordinates'] = self.coordinates
        dump['variables'] = dict(
            longitude = self.variables['longitude'].json_dict(),
            latitude = self.variables['latitude'].json_dict(),
            time = self.variables['time'].json_dict(),
            height = self.variables['height'].json_dict()
        )
        for item in var:
            dump['variables'][item] = self.variables[item].json_dict()
        url = '%s/%s' % (self.params['path'], '.atlantis')
        f = open(url, 'w')
        json.dump(dump, f, indent=2)
        f.close()

        self.__init__(path=self.params['path'])

        return


    def refresh(self, var=[]):
        """Refreshes the Atlantis dataset."""
        raise ValueError('This function is not implemented yet.')
        return


    def read(self, t=None, z=None, y=None, x=None, N=None, K=None, J=None,
        I=None, var=None, nonan=True, result='full', merge=False,
        profile=False, dummy=False):
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
            merge (boolean, optional) :
                If true, and multiple variables are selected, then merges
                the data by summation. Default is false.
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
        # Switching to raw dataset indices
        II, JJ = self.params['lon_i'][JJ, II], self.params['lat_j'][JJ, II]
        # Initializes return variable
        if var == None:
            var_list = list(set(self.variables.keys()) -
                set(self.default_vars))
        elif type(var) == str:
            var_list = [var]
        else:
            var_list = list(var)
        var = dict()
        if merge:
            var['merged'] = ma.zeros(shape)
        else:
            for nvar in var_list:
                var[nvar] = ma.zeros(shape)
        var_list_length = len(var)
        if DEBUG:
            print '\r\n',
            print 't: ', t
            print 'z: ', z
            print 'y:', y
            print 'x:', x
            print 'N: ', N
            print 'K: ', K
            print 'J: ', J
            print 'I:', I
            print 'shape: ', shape
            print 'var_list', var_list
        # Walks through every time index and loads data range from maps.
        for n, T in enumerate(t):
            t2 = time()
            if profile:
                s = '\rLoading \'{0}\' data... {1} '.format(nvar,
                    profiler(shape[0], n + 1, 0, t1, t2))
                self._message(s)
            for nvar in var_list:
                # Reads numpy map file, masks where NaN and if appropriate,
                # change masked data values to zero.
                fname = '%s/%s' % (self.params['path'],
                    self.params['file_list_{}'.format(nvar)][N[n]])
                #try:
                if True:
                    P = loadtxt(fname)[JJ+1, II+1]
                #except:
                else:
                    print fname, shape
                    P = ones(shape[2:]) * nan
                P = ma.masked_where(isnan(P), P)
                try:
                    if nonan and (type(P.mask) != bool_):
                        P.data[P.mask] = 0
                except:
                    pass
                #
                if merge:
                    var['merged'][n, 0, :, :] = (var['merged'][n, 0, :, :] +
                        P[None, None, :, :])
                else:
                    var[nvar][n, 0, :, :] = P[None, None, :, :]
            #
            if DEBUG:
                if (n+1) % 10 == 0:
                    self._message('|')
                elif (n+1) % 5 == 0:
                    self._message(':')
                else:
                    self._message('.')

        if profile:
            self._message('\r\n')

        if var_list_length == 1:
            var = var[var.keys()[0]]

        if DEBUG:
            print
            print 'var: ', var

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


    def write(self, t, z, var_name=None, prefix='', suffix='', fmt='%.18e'):
        """Writes the data."""
        # Checks for dimension of input data
        if z.shape != (self.dimensions['j'], self.dimensions['i']):
            raise ValueError('Input data dimensions do not match dataset '
                'dimensions.')
        # If no variable name is set, then uses first variable from variable
        # list
        var_list = list(set(self.variables.keys()) - set(self.default_vars))
        n_var = len(var_list)
        if var_name == None:
            if n_var == 0:
                var_name = var_list[0]
            else:
                raise ValueError('For multivariate datasets, variable name '
                    'has to be specified.')

        if n_var > 1:
            fpath = '{0}/{1}'.format(self.params['path'], var_name)
        else:
            fpath = self.params['path']

        # Converts masked array to numpy array making sure that masked values
        # become NaN's
        if type(z) == ma.MaskedArray:
            z.data[z.mask] = nan
            z = z.data

        dat = zeros((self.dimensions['j'] + 1, self.dimensions['i'] + 1))
        dat[0, 0] = t
        dat[0, 1:] = self.variables['longitude'].data
        dat[1:, 0] = self.variables['latitude'].data
        dat[1:, 1:] = z
        url = '%s/%s%s_%06d%s.xy.gz' % (fpath, prefix, var_name, t, suffix)
        savetxt('%s' % (url), dat, fmt=fmt, delimiter='\t')

        return


    def climatology(self, var=None, w=None, result='year', profile=True,
        **kwargs):
        """Returns monthly climatology of the time-series.

        PARAMETERS
            var (string, optional) :
                Indicates which variable of the grid will be read. If
                the parameter is a list of variables, then the data
                will be returned as a list of arrays.
            w (array like) :
                Data weight, should have same dimensions as 'z'.
            result (string) :
                If set to 'year', returns only one year of data. If set
                to 'full', returns the climatology for every time t.

        OPTIONAL KEYWORD ARGUMENTS
            z, y, x (array like, optional) :
                Sets the height, latitude and longitude for which the
                data will be read.
            K, J, I (array like, optional) :
                Sets the vertical, meridional and zonal indices for
                which the data will be read.

        RETURN
            var_clim (array like) :
                Climatological averages.

        """
        t1 = time()
        # Checks if optional keywords are fine
        exargs = set(kwargs) - set(['z', 'y', 'x', 'K', 'J', 'I'])
        if exargs:
            raise Warning('Invalid arguments: %s' % list(exargs))
        kwargs['result'] = 'var only'

        # Checks for proper dimensions loading dummy data.
        d, c, b, a = self.read(dummy=True, **kwargs)
        if w != None:
            if w.shape != (d, c, b, a):
                raise Warning, 'Data and weight arrays are not the same.'

        # Get list of variables to be loaded.
        if var == None:
            var = self.params['var_list']

        # Starts converting time to datetime format. Determines the start and end
        # of the relevant dataset to ensure that only whole years are used.
        # Initializes climatology variable and calculates the averages.
        start, end, Time = self.getFullYearsIndices()
        time_clim = arange(12) + 1
        var_clim = dict()
        for i in range(12):
            t2 = time()
            if profile:
                s = '\rCalculating climatology... %s ' % (profiler(12, i, 0,
                    t1, t2),)
                self._message(s)
            # Initializes sum and weight arrays, selects appropriate data range
            # and walksthrough every data map (grid) to calculate
            # climatological mean.
            var_sum, var_weight = dict(), dict()
            selt = flatnonzero(Time[start:end+1, 1] == i + 1) + start
            n = len(selt)
            for j, jj in enumerate(selt):
                z = self.read(N=[jj], var=var, **kwargs)
                if type(z) != dict:
                    z = dict(var=z)
                # Walks through each item in z dictionary. This is wise, since
                # the grid class may return keys different from the variable
                # list 'var', i.e. 'tauxy' if 'taux' and 'tauy' are selected.
                for item in z.keys():
                    if item not in var_clim.keys():
                        var_clim[item] = ma.zeros([12, c, b, a],
                            dtype=z[item].dtype)
                    if item not in var_sum.keys():
                        var_sum[item] = ma.zeros([c, b, a],
                            dtype=z[item].dtype)
                        var_weight[item] = ma.zeros([c, b, a],
                            dtype=z[item].dtype)
                    #
                    w = ~z[item][0, :, :, :].mask + 0.
                    var_sum[item] += z[item][0, :, :, :].data * w
                    var_weight[item] += w
                #
                if profile:
                    s = '\rCalculating climatology... %s ' % (profiler(12,
                        i + (j + 1) / n,  0,  t1, t2),)
                    self._message(s)
            #
            for item in var_clim.keys():
                var_clim[item][i, :, :, :] = var_sum[item] / var_weight[item]
        #
        if profile:
            self._message('\r\n')
        #
        for item in var_clim.keys():
            var_clim[item].mask = (var_clim[item].mask |
                isnan(var_clim[item].data))
            if result == 'full':
                var_clim[item] = var_clim[item][Time[:, 1]-1]
        #
        if len(var_clim.keys()) == 1:
            var_clim = var_clim[var_clim.keys()[0]]
        #
        start = '%04d-%02d' % (Time[start][0], Time[start][1])
        end = '%04d-%02d' % (Time[end][0], Time[end][1])
        #
        return var_clim, Time, start, end


    def getFullYearsIndices(self, offset=None):
        """
        Returns the time indices for beginning and end of full year
        time-series. In case the dataset has no time attribute, then
        it will consider one year of data.

        PARAMETER
            offset (float, optional) :
                Sets the offset for full years.

        RETURNS
            start, end (integer) :
                Indexes of time-series.
            Time (array) :
                Array containing year, month, day, hour and minute.

        """
        if offset != None:
            raise Warning('Option not implemented yet.')

        if self.variables['time'].data == None:
            return 0, 364

        # Starts converting time to datetime format. Determines the start and end
        # of the relevant dataset to ensure that only whole years are used.
        Time = num2ymd(self.variables['time'].data)
        try:
            start = flatnonzero(Time[:, 1] == 1)[0]
        except:
            start = None
        try:
            end = flatnonzero(Time[:, 1] == 12)[-1]
        except:
            end = None

        return start, end, Time


    def getLongitudeLatitudeLimits(self, lon, lat, xlim=None, ylim=None):
        """
        Returns the spatial indices according to the longitude
        and latitude limits.

        Parameters
        ----------
        lon, lat : array like
            Longitude and latitude arrays.
        xlim, ylim : list, optional
            Longitude and latitude limits.

        Returns
        -------
        lon, lat : array like
            Longitude and latitude arrays according to limits.
        xlim, ylim : array like
            Sorted limits.
        ii, jj : array like
            Mesh arrays of zonal and meridional indices.

        """
        # If xlim and ylim are set, calculate how many indices have to be moved
        # in order for latitude array to start at xlim[0].
        if xlim == None:
            xlim = (lon.min(), lon.max())
        if ylim == None:
            ylim = (lat.min(), lat.max())
        #
        LON = lon_n(lon, xlim[1])
        i = argsort(LON)
        selx = i[flatnonzero((LON[i] >= xlim[0]) & (LON[i] <= xlim[1]))]
        sely = flatnonzero((lat >= ylim[0]) & (lat <= ylim[1]))
        ii, jj = meshgrid(selx, sely)
        lon = LON[selx]
        lat = lat[sely]
        #
        return lon, lat, xlim, ylim, ii, jj


    def createMeterGrid(self, units='km'):
        """Creates a meter grid based on longitude and latitude degrees."""
        self.variables['xm'] = Variable(canonical_units=units,
            description='Zonal distance.')
        self.variables['ym'] = Variable(canonical_units=units,
            description='Meridional distance.')
        self.variables['xm'].data, self.variables['ym'].data = (
            metergrid(self.variables['longitude'].data,
            self.variables['latitude'].data, units=units)
        )


    def derivative(self, A, axis='longitude', p=1, q=3, cyclic=True,
        units=None, mask=None):
        """Higher order differential.

        Calculates the p-th derivative of `A` with respect to the given
        axis using stencils of width n. As default, it assumes that `A` is a
        mapped field on the globe, so it is periodic in the x-direction.
        It exploits this peridicity so that there are not missing data
        points at the boundaries.

        The gradient is computed using central differences in the
        interior and first differences at the boundaries. The returned
        gradient hence has the same shape as the input array.

        PARAMETERS
            A (array like) :
                Data to calculate the derivate. It has to have the same
                dimensions as the data grid.
            axis (string, optional) :
                Axis onto which the derivative will be calculated.
                Values should be either 'longitude' (default), 'xm',
                'latitude', 'ym'.
            p (integer, optional) :
                Order of the derivative to be calculated. Default is to
                calculate the first derivative (p=1). 2*n-p+1 gives the
                relative order of the approximation.
            q (integer, optional) :
                Length of the stencil used for centered differentials.
                The length has to be odd numbered. Default is q=3.
            cyclic (boolean, optional) :
                Sets whether `A` is to be considered zonally periodic.
                Default is true.
            units (string, optional) :
                Unit to which the spatial coordinates will be converted.
                Accepted values are 'm' for meters, 'km' for kilometers,
                and 'nm' for International nautical miles or 'deg' for
                degrees. If not set (None), same units as axis will be
                used.
            mask (array like, optional) :

        RETURNS
            dA (array like) :
                The calculated derivate with same dimensions as A

        REFERENCES
            Cushman-Roisin, B. & Beckers, J.-M. Introduction to
            geophysical fluid dynamics: Physical and numerical aspects.
            Academic Press, 2011, 101, 828.

            Arbic, Brian B. Scott, R. B.; Chelton, D. B.; Richman, J. G.
            and Shriver, J. F. Effects of stencil width on surface ocean
            geostrophic velocity and vorticity estimation from gridded
            satellite altimeter data. Journal of Geophysical Research,
            2012, 117, C03029.

        """
        global DEBUG
        # Checks the dimensions of the input and mask arrays
        if A.shape != (self.dimensions['j'], self.dimensions['i']):
            raise ValueError('Input data dimensions do not match dataset '
                'dimensions.')
        if mask != None:
            if mask.shape != (self.dimensions['j'], self.dimensions['i']):
                raise ValueError('Mask arraydimensions do not match dataset '
                    'dimensions.')

        # Checks spacial units
        ax_units = self.variables[axis].canonical_units
        if ax_units in ['degree_north', 'degree_south', 'degree_east',
            'degree_west']:
            ax_units = 'deg'
        scale = 1.
        if (units != None) & (ax_units != units):
            if (ax_units == 'km') & (units == 'm'):
                scale = 1e-3
            elif (ax_units == 'm') & (units == 'km'):
                scale = 1e3
            elif (ax_units == 'deg') & (units in ['m', 'km', 'nm']):
                # Checks if meter grid exists in data set. If not,
                # generates it.
                if len(set(['xm', 'ym']) - set(self.variables.keys())) > 0:
                    self.createMeterGrid(units=units)
                    if axis == 'longitude':
                        axis = 'xm'
                    elif axis == 'latitude':
                        axis = 'ym'
            else:
                raise ValueError('Unit conversion from \'%s\' to \'%s\' not '
                    'implemented yet.' % (ax_units, units))

        # Determines the axis direction
        if axis in ['longitude', 'xm']:
            direction = 'i'
        elif axis in ['latitude', 'ym']:
            direction = 'j'
        else:
            raise ValueError('Invalid axis \'%s\'.' % (axis))

        # Makes shure the length of the stencil is odd numbered.
        q += (q % 2) - 1
        # Calculate left and right stencils.
        q_left = (q - 1) / 2
        q_right = (q - 1) / 2

        # If data is cyclic, pad data on boundaries according to the stencile
        # width.
        if cyclic:
            A = hstack([A[:, -q_left:], A, A[:, :q_right]])
            suffix = '_c'
        else:
            suffix = ''
        J, I = A.shape

        # Stencil coefficients unique ID.
        coeffs_id = '%s_p%d_q%d%s' % (axis, p, q, suffix)

        # If any parameter has changed or if the stencil coefficients for the
        # chosen direction has not been calculated yet, then calculate them.
        if coeffs_id not in self.stencil_coeffs.keys():
            # Set axis variable and do some speed optimization if values are
            # repeated along direction axis.
            axis = self.variables[axis].data
            if axis.ndim == 1:
                if direction == 'i':
                    axis = axis[None, :]
                elif direction == 'j':
                    axis = axis[:, None]
            if cyclic & (axis.shape[1] > 1):
                axis = hstack([(axis[:, 0, None] - axis[:, -1, None]) *
                    ones([1, q_left]) + axis[:, -q_left-1:-1], axis,
                    (axis[:, -1, None] - axis[:, 0, None]) *
                    ones([1, q_right]) + axis[:, 1:q_right+1]])
            #
            self.stencil_coeffs[coeffs_id] = (
                self._derivative_stencil_coefficients(axis, direction, p=p,
                    q=q, cyclic=cyclic)
            )

        # Calculate the p-th derivative
        dA = zeros((J, I))
        coeffs = self.stencil_coeffs[coeffs_id]
        for i in arange(q) - q_left:
            if direction == 'i':
                if i < 0:
                    u, v = -i, I
                else:
                    u, v = 0, I - i
                dA[:, u:v] += (coeffs[:, u:v, i+q_left] * A[:, u+i:v+i])
            elif direction == 'j':
                if i < 0:
                    u, v = -i, J
                else:
                    u, v = 0, J - i
                dA[u:v, :] += (coeffs[u:v, :, i+q_left] * A[u+i:v+i, :])

        if cyclic:
            return dA[:, q_left:-q_right] * scale**p
        else:
            return dA * scale**p


    def _derivative_stencil_coefficients(self, axis, direction, p=1, q=3,
        cyclic=True, mask=None, profile=True):
        """Calculates the coefficients needed for the derivative.

        PARAMETERS
            axis (string) :
                Axis onto which the derivative will be calculated.
                Values should be either 'longitude' (default), 'xm',
                'latitude', 'ym'.
            direction (char) :
            p (integer, optional) :
                Order of the derivative to be calculated. Default is to
                calculate the first derivative (p=1). 2*n-p+1 gives the
                relative order of the approximation.
            q (integer, optional) :
                Length of the stencil used for centered differentials.
                The length has to be odd numbered. Default is q=3.
            cyclic (boolean, optional) :
            mask (array like, optional) :
            profile (boolean, optional) :

        RETURNS
            Coefficients (a_q) needed for the linear combination of `q`
            points to get the first derivative according to Arbic et
            al. (2012) equations (20) and (22). At the boundaries
            forward and backward differences approximations are
            calculated.

        REFERENCES
            Cushman-Roisin, B. & Beckers, J.-M. Introduction to
            geophysical fluid dynamics: Physical and numerical aspects
            Academic Press, 2011, 101, 828

            Arbic, Brian B. Scott, R. B.; Chelton, D. B.; Richman, J.
            G. & Shriver, J. F. Effects of stencil width on surface
            ocean geostrophic velocity and vorticity estimation from
            gridded satellite altimeter data. Journal of Geophysical
            Research, 2012, 117, C03029

        """
        global DEBUG
        t1 = time()
        #
        descriptor = 'p=%d, q=%d, direction=%s, cyclic=%s' % (p, q, direction,
            cyclic)
        # Calculate left and right stencils.
        q_left = (q - 1) / 2
        q_right = (q - 1) / 2 + 1
        #
        J, I = axis.shape
        #
        if mask == None:
            mask = ones([J, I])
        else:
            mask = 1. * array(~mask)

        # Constructs matrices according to Cushman-Roisin & Beckers (2011)
        # equations (1.25) and adapted for variable grids as in Arbic et
        # al. (2012), equations (20), (22). The linear system of equations
        # is solved afterwards.
        coeffs = zeros((J, I, q))
        smart_coeffs = dict()
        for j in range(J):
            if profile:
                t2 = time()
                s = '\rCalculating coefficients, %s... %s ' % (descriptor,
                    profiler(J, j, 0, t1, t2))
                self._message(s)
            for i in range(I):
                A = zeros((q, q))
                if direction == 'i':
                    if i < q_left:
                        start = q_left - i
                    else:
                        start = 0
                    if i > I - q_right:
                        stop = i - (I - q_right)
                    else:
                        stop = 0
                elif direction == 'j':
                    if j < q_left:
                        start = q_left - j
                    else:
                        start = 0
                    if j > J - q_right:
                        stop = j - (J - q_right)
                    else:
                        stop = 0
                A[0, start:q+stop] = 1
                if direction == 'i':
                    da = axis[j, i-q_left+start:i+q_right-stop] - axis[j, i]
                elif direction == 'j':
                    da = axis[j-q_left+start:j+q_right-stop, i] - axis[j, i]
                da_key = str(da)
                if da_key not in smart_coeffs.keys():
                    for h in range(1, q):
                        A[h, start:q-stop] = da ** h
                    B = zeros((q, 1))
                    # This tells where the p-th derivative is calculated
                    B[p] = factorial(p)
                    C = linalg.solve(A[:q-start-stop, start:q-stop],
                        B[:q-(start+stop), :])
                    smart_coeffs[da_key] = C.flatten()
                    if DEBUG:
                        print 'axis =', axis
                        print 'j = %d, i = %d:' % (j, i)
                        print 'A = ', A
                        print 'B = ', B
                        print 'C = ', C
                coeffs[j, i, start:q-stop] = smart_coeffs[da_key]
            #
            if profile:
                s = '\rCalculating coefficients, %s... %s ' % (descriptor,
                    profiler(J, j+1, 0, t1, t2))
                self._message(s)
        #
        if profile:
            self._message('\n\r')
        #
        if DEBUG:
            print 'coeffs[%s%d] = ' % (direction, p), coeffs, coeffs.shape

        # Reshapes the coefficient if the dimension is one
        if direction == 'i':
            if J == 1:
                coeffs = coeffs.repeat(self.dimensions['j'], axis=0)
        elif direction == 'j':
            if I == 1:
                if cyclic:
                    i0 = q - 1
                else:
                    i0 = 0
                coeffs = coeffs.repeat(self.dimensions['i'] + i0, axis=1)

        return coeffs


    def derivative_smoother(self, A, q=3, n=15, epsilon=1e-5):
        """Smoothes a scalar array through derivatives and integrals.

        PARAMETERS
            A (array like) :
            q (integer, optional) :
                Length of the stencil used for centered differentials.
                The length has to be odd numbered. Default is q=3.
            n (integer, optional) :
                Number of maximum iterations.
            epsilon (float, optional) :
                Root mean square difference between subsequent
                iterations.

        RETURNS
            Smoothed scalar array.

        """
        global DEBUG
        # Checks the dimensions of the input array
        if A.shape != (self.dimensions['j'], self.dimensions['i']):
            raise ValueError('Input data dimensions do not match dataset '
                'dimensions.')

        # Checks if input array is of masked type
        if type(A) == ma.MaskedArray:
            masked = True
            mask = A.mask
            A = A.data
        else:
            masked = False

        # Determines zonal and meridional grid spacing
        try:
            dx = self.params['dlon']
            dy = self.params['dlat']
        except:
            dx = (self.variables['longitude'].data[1] -
                self.variables['longitude'].data[0])
            dy = (self.variables['latitude'].data[1] -
                self.variables['latitude'].data[0])

        # Keeps calculating the first derivative of A and integrate it back
        # until maximum number of iterations 'n' are executed or until RMSD is
        # less than 'epsilon'.
        for i in range(n):
            if i == 0:
                a = A
            else:
                a = smooth_A
            dadx = self.derivative(a, axis='longitude', p=1, q=q)
            a_left = a[:-1, :] + dadx[1:, :] * dx
            a_right = a[1:, :] - dadx[:-1, :] * dx
            dady = self.derivative(a, axis='latitude', p=1, q=q)
            a_top = a[:, :-1] + dady[:, 1:] * 0.25
            a_bottom = a[:, 1:] - dady[:, :-1] * 0.25
            #
            w_sum = zeros(A.shape)
            a_sum = zeros(A.shape, dtype=float)
            a_sum[1:, :] += a_left; w_sum[1:, :] += 1
            a_sum[:-1, :] += a_right; w_sum[:-1, :] += 1
            a_sum[:, 1:] += a_top; w_sum[:, 1:] += 1
            a_sum[:, :-1] += a_bottom; w_sum[:, :-1] += 1
            smooth_A = a_sum / w_sum
            smooth_A = ma.masked_where(isnan(smooth_A), smooth_A)
            #
            try:
                delta = ((a - smooth_A)**2.).mean() - RMSD
            except:
                delta = -65535.
            RMSD = ((a - smooth_A)**2.).mean()
            if DEBUG:
                print '%d. RMSD=%f, delta=%f' % (i, RMSD, delta), (RMSD <= epsilon)
            if RMSD <= epsilon:
                break
            if delta > 0.:
                smooth_A = a
                break

        if masked:
            smooth_A = ma.masked_where(mask, smooth_A)
        #
        return smooth_A


    def curl(self, A, axis_x='longitude', axis_y='latitude', **kwargs):
        """Returns the curl of a two-dimensional vector field according
        to a given stencil width.

        The input vector array should be in complex notation, such that
        A = u + j * v, where j = (-1)**0.5.

        Calculates the first derivative of `A` with respect x and y using
        q stencils. As default, it assumes that `A` is a mapped vector
        field on the globe, so it is periodic in the x-direction. It
        exploits this peridicity so that there are no missing data points
        at the boundaries.

        The curl is computed using central differences in the interior
        and first differences at the boundaries. The returned gradient
        hence has the same shape as the input data array.

        PARAMETERS
            A (array like) :
                Complex two-dimensional data array.
            axis_x, axis_y (string, optional) :
                The two axis onto which the derivative will be
                calculated. Values should be either 'longitude'
                (default), 'xm', 'latitude' (default), 'ym'.
            q (integer, optional) :
                Length of the stencil used for centered differentials.
                The length has to be odd numbered. Default is n=3.
            cyclic (boolean, optional) :
                Sets whether `A` is to be considered zonally periodic.
                Default is true.
            units (string, optional) :
                Unit to which the spatial coordinates will be converted.
                Accepted values are 'm' for meters, 'km' for kilometers,
                and 'nm' for International nautical miles or 'deg' for
                degrees. If not set (None), save units as axis will be
                used.

        RETURNS
            C (ndarray) :
                Arrays of the same shape as `A` giving the curl of `A`.

        """
        dvdx = self.derivative(A.imag, axis=axis_x, p=1, **kwargs)
        dudy = self.derivative(A.real, axis=axis_y, p=1, **kwargs)
        return dvdx - dudy


    def grad(self, A, axis_x='longitude', axis_y='latitude', **kwargs):
        """Calculates the horizontal gradient of a two-dimensional
        scalar field.

        \\nabla A = \frac{\partial A}{\partial x}\hat{\imath} +
            \frac{\partial A}{\partial y} \hat{\jmath}

        PARAMETERS
            A (array like) :
                Complex two-dimensional data array.
            axis_x, axis_y (string, optional) :
                The two axis onto which the derivative will be
                calculated. Values should be either 'longitude'
                (default), 'xm', 'latitude' (default), 'ym'.
            q (integer, optional) :
                Length of the stencil used for centered differentials.
                The length has to be odd numbered. Default is n=3.
            cyclic (boolean, optional) :
                Sets whether `A` is to be considered zonally periodic.
                Default is true.
            units (string, optional) :
                Unit to which the spatial coordinates will be converted.
                Accepted values are 'm' for meters, 'km' for kilometers,
                and 'nm' for International nautical miles or 'deg' for
                degrees. If not set (None), save units as axis will be
                used.

        RETURNS
            C (ndarray) :
                Arrays of the same shape as `A` giving the gradient of A
                in complex notation, such that the zonal component is
                the real part and the meridional component is the
                imaginary part.

        """
        C = (self.derivative(A, axis=axis_x, p=1, **kwargs) +
            1j * self.derivative(A, axis=axis_y, p=1, **kwargs))
        return C


    def div(self, A, axis_x='longitude', axis_y='latitude', **kwargs):
        """Calculates the divergence of a two-dimensional vector field.

        The input vector array should be in complex notation, such that
        A = u + j * v, where j = (-1)**0.5.

        PARAMETERS
            A (array like) :
                Complex two-dimensional data array.
            axis_x, axis_y (string, optional) :
                The two axis onto which the derivative will be
                calculated. Values should be either 'longitude'
                (default), 'xm', 'latitude' (default), 'ym'.
            q (integer, optional) :
                Length of the stencil used for centered differentials.
                The length has to be odd numbered. Default is n=3.
            cyclic (boolean, optional) :
                Sets whether `A` is to be considered zonally periodic.
                Default is true.
            units (string, optional) :
                Unit to which the spatial coordinates will be converted.
                Accepted values are 'm' for meters, 'km' for kilometers,
                and 'nm' for International nautical miles or 'deg' for
                degrees. If not set (None), save units as axis will be
                used.

        RETURNS
            C (ndarray) :
                Arrays of the same shape as `A` giving the divergence
                of `A`.

        """
        C = (self.derivative(A.real, axis=axis_x, p=1, **kwargs) +
            self.derivative(A.imag, axis=axis_y, p=1, **kwargs))
        return C


    def laplacian(self, A, axis_x='longitude', axis_y='latitude', **kwargs):
        """Calculates the Laplacian of a two-dimensional scalar field.

        The Laplace operator in two dimesions is given by

            \\nabla^2 = \frac{\partial^2 A}{\partial x^2} +
                \frac{\partial^2 A}{\partial y^2}

        PARAMETERS
            A (array like) :
                Complex two-dimensional data array.
            axis_x, axis_y (string, optional) :
                The two axis onto which the derivative will be
                calculated. Values should be either 'longitude'
                (default), 'xm', 'latitude' (default), 'ym'.
            q (integer, optional) :
                Length of the stencil used for centered differentials.
                The length has to be odd numbered. Default is n=3.
            cyclic (boolean, optional) :
                Sets whether `A` is to be considered zonally periodic.
                Default is true.
            units (string, optional) :
                Unit to which the spatial coordinates will be converted.
                Accepted values are 'm' for meters, 'km' for kilometers,
                and 'nm' for International nautical miles or 'deg' for
                degrees. If not set (None), save units as axis will be
                used.

        RETURNS
            C (ndarray) :
                Arrays of the same shape as `A` giving the divergence
                of `A`.

        TODO
            Convert the Laplacian operator on a sphere!

        """
        C = (self.derivative(A, axis=axis_x, p=2, **kwargs) +
            self.derivative(A, axis=axis_y, p=2, **kwargs))
        return C


    @property
    def name(self):
        """Returns a human-readable name."""
        return self._name

    @name.setter
    def name(self, s='Unknown'):
        self._name = s

    @property
    def description(self):
        return self._description

    @description.setter
    def description(self, s):
        self._description = s


def get_standard_variable(name, **kwargs):
    """
    Reads the CF standard name table xml file and returns the default
    parameters for the standard variable set in 'name'. It the desired
    variable is an alias, standard name is retured.

    """
    global DEBUG

    # 1. Creates a variable instance and sets its attributes.
    var = Variable()
    var.standard_name = None
    var.canonical_units = None
    var.description = None
    var.grib = None
    var.amip = None
    # 2. Locates and opens the cf-standard-name-table.xml file
    _path = path.dirname(__file__)
    _fname = 'cf-standard-name-table.xml'
    xml = etree.parse('%s/%s' % (_path, _fname))
    # 3. Reads the xml document, searches for chosen standard variable and
    # sets variable attributes. First looks at aliases and then at standard
    # names.
    root = xml.getroot()
    try:
        entry = root.xpath("//alias[@id='%s']" % (name))[0]
        for item in entry.getchildren():
            if item.tag == 'entry_id':
                name = item.text
    except:
        pass
    try:
        entry = root.xpath("//entry[@id='%s']" % (name))[0]
        var.standard_name = entry.attrib['id']
        for item in entry.getchildren():
            if DEBUG:
                print 'get_standard_variable: %s = %s' % (item.tag, item.text)
            if item.tag == 'canonical_units':
                var.canonical_units = item.text
            elif item.tag == 'description':
                var.description = item.text
            elif item.tag == 'grib':
                var.grib = item.text
            elif item.tag == 'amip':
                var.amip = item.text
    except:
        raise ValueError("'%s' is not a valid standard name." % (name))
    # 5. Sets attribute from input arguments
    for key, item in kwargs.items():
        setattr(var, key, item)
    # 4. Returns standard variable
    return var
