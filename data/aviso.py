# -*- coding: utf-8 -*-
"""Atlantis data framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional data and
associated metadata. Much of this is based uppon Iris library

This module implements Aviso merged sea level anomaly reference series
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
    1 (2013-11-29 16:00 -0300 DST)

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

from matplotlib import dates
from numpy import (arange, argsort, array, asarray, flatnonzero, in1d, isnan, 
    ma, meshgrid, nan, ones)
from os import listdir, remove
from os.path import isfile
from sys import stdout
from time import time
from gzip import open as gzopen
from uuid import uuid1 as uuid
from glob import glob
from scipy.io import netcdf_file as netcdf

import atlantis.data

from klib.common import profiler, lon_n, reglist
from atlantis.astronomy import metergrid

DEBUG = False

class Grid(atlantis.data.Grid):
    """Common grid for Aviso merged gridded delayed time products (sea
    level anomaly, geostrophic currents and errors).
    
    """
    def __init__(self, path=None, mask_file=None, xlim=None, ylim=None,
        tlim=None, useqd=False):
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
        self.stencil_coeffs = dict()
        self.stencil_params = dict()
        if useqd:
            self.params['datasets'] = [dict(id='h', var='h_qd')]
            self.params['var_dict'] = dict(h_Grid_0001 = 'h_qd')
            self.params['var_tcid'] = dict(h_qd=['h', 'h_qd', 'Grid_0001'])
        else:
            self.params['datasets'] = [dict(id='h', var='h'),
                dict(id='uv', var='uv'), dict(id='err', var='err')]
            self.params['var_dict'] = dict(
                h_Grid_0001 = 'h',
                uv_Grid_0001 = 'u',
                uv_Grid_0002 = 'v',
                err_Grid_0001 = 'err'
            )
            self.params['var_tcid'] = dict(
                h = ['h', 'h', 'Grid_0001'],
                u = ['uv', 'uv', 'Grid_0001'],
                v = ['uv', 'uv', 'Grid_0002'],
                err = ['err', 'err', 'Grid_0001']
            )
        # Creates an universally unique identifiers (UUID) for this instance
        self.params['uuid'] = str(uuid())
        
        # Sets global parameters for grid.
        if path == None:
            path = ('/home/sebastian/academia/data/aviso/msla/merged')
        self.params['path'] = path
        self.params['mask_file'] = mask_file
        self.params['missing_value'] = -9999.
        
        # Generates list of files, tries to match them to the pattern and to 
        # extract the time.
        file_pattern = ('dt_ref_global_merged_msla_(%s)_(\d*)_(\d*)_(\d*)'
            '.nc.gz' % ('|'.join([item['var'] for item in
            self.params['datasets']])))
        flist = listdir('%s/%s' % (self.params['path'],
            self.params['datasets'][0]['id']))
        flist.sort()
        flist, match = reglist(flist, file_pattern)
        
        # Convert dates to matplotlib format, i.e. days since 0001-01-01 UTC.
        time_list = array(dates.datestr2num(['%4s-%2s-%2s 12:00' %
            (item[1][:4], item[1][4:6], item[1][6:]) for item in match]))
        
        # If tlim are set, calculate the time limits of the dataset and
        # corresponding files.
        if tlim != None:
            for i, t in enumerate(tlim):
                if type(t) == str:
                    tlim[i] = dates.datestr2num(t)
            #
            t_sel = flatnonzero(((time_list >= tlim[0]) &
                (time_list <= tlim[1])))
            time_list = time_list[t_sel]
        else:
            t_sel = range(len(time_list))
        
        fdict = [dict(start=match[n][1], end=match[n][2],
            creation=match[n][3]) for n in t_sel]
        self.params['file_list'] = fdict
        if len(flist) == 0:
            return
        
        # Reads first file in dataset to determine array geometry and 
        # dimenstions (lon, lat)
        params = dict(path=self.params['path'],
            dataset=self.params['datasets'][0]['id'],
            datavar=self.params['datasets'][0]['var'],
            **self.params['file_list'][0])
        fname = self.create_filename(**params)
        data = self.read_file(fname)
        lat = data.variables['NbLatitudes'].data
        lon = data.variables['NbLongitudes'].data
        
        # If xlim and ylim are set, calculate how many indices have to be moved
        # in order for latitude array to start at xlim[0].
        lon, lat, xlim, ylim, ii, jj = self.getLongitudeLatitudeLimits(lon,
            lat, xlim, ylim)
        self.params['xlim'], self.params['ylim'] = xlim, ylim
        self.params['lon_i'], self.params['lat_j'] = ii, jj
        self.params['dlon'] = lon[1] - lon[0]
        self.params['dlat'] = lat[1] - lat[0]
        
        # Initializes the grid attributes, dimensions, coordinates and
        # variables.
        self.name = 'sea_level_anomaly_geostrophic_velocities'
        for attr, attr_value in vars(data).iteritems():
            if attr in ['mode', 'filename']:
                continue
            if type(attr_value) == str:
                if attr in ['name']:
                    self.name = attr_value
                elif attr in ['description', 'summary', 'title']:
                    self.description = attr_value
                else:
                    self.attributes[attr.lower()] = attr_value
        self.dimensions = dict(n=time_list.size, k=1, j=lat.size, i=lon.size)
        self.coordinates = dict(n='time', k='height', j='latitude',
            i='longitude')
        #
        self.variables = dict(
            time = atlantis.data.variable(
                canonical_units='days since 0001-01-01 UTC',
                data=time_list,
            ),
            height = atlantis.data.get_standard_variable('height', data=[0.]),
            latitude = atlantis.data.get_standard_variable('latitude',
                data=lat),
            longitude = atlantis.data.get_standard_variable('longitude',
                data=lon),
            xm = atlantis.data.variable(
                canonical_units = 'km',
                description = 'Zonal distance.'
            ),
            ym = atlantis.data.variable(
                canonical_units = 'km',
                description = 'Meridional distance.'
            ),
        )
        #
        self.variables['xm'].data, self.variables['ym'].data = (
            metergrid(self.variables['longitude'].data, 
            self.variables['latitude'].data, unit='km')
        )
        # Walks through every dataset to read list of variables.
        self.params['var_list'] = list()
        for i, dataset in enumerate(self.params['datasets']):
            if i > 0:
                params = dict(path=self.params['path'],
                    dataset=dataset['id'], datavar=dataset['var'],
                    **self.params['file_list'][0])
                fname = self.create_filename(**params)
                data = self.read_file(fname)
            # Walks through every variable in NetCDF file
            for var in data.variables.keys():
                if var in ['Grid_0001', 'Grid_0002']:
                    nvar = self.params['var_dict']['{0}_{1}'.format(
                        dataset['id'], var)]
                    attribs = dict(
                        missing_value = data.variables[var]._FillValue,
                        canonical_units = data.variables[var].units,
                        description = data.variables[var].long_name,
                        dataset = dataset,
                        variable = var
                    )
                    self.variables[nvar] = atlantis.data.variable(**attribs)
                    self.params['var_list'].append(nvar)
            # Closes the data access and removes temporary NetCDF file
            self.close_file(data)
        
        return
    
    
    def create(*args, **kwargs):
        return
    
    
    def create_filename(*args, **kwargs):
        """Determines the file name according to the input parameters.

        PARAMETERS
            path (string, optional) :
            dataset (string, optional) :
            datavar (string, optional) :
            start (string, optional) :
            end (string, optional) :
            creation (string, optional) :

        RETURNS
            filename (string) :

        """
        default_args = dict(path='', dataset='', datavar='', start='', end='',
            creation='')
        kwargs_keys = kwargs.keys()
        for key, item in default_args.items():
            if key not in kwargs_keys:
                kwargs[key] = item
        fname = ('{path}/{dataset}/dt_ref_global_merged_msla_{datavar}_'
            '{start}_{end}_{creation}.nc.gz').format(**kwargs)

        if isfile(fname):
            return fname
        else:
            pattern = ('{path}/{dataset}/dt_ref_global_merged_msla_{datavar}_'
            '{start}_{end}_*.nc.gz').format(**kwargs)
            flist = glob(pattern)
            if len(flist) > 0:
                return flist[-1]
            else:
                raise ValueError('No matching file found.')
    
    
    def read_file(self, filename):
        """Reads zipped NetCDF file and returns its file pointer.
        
        """
        # Uncompress NetCDF file.
        f = gzopen('%s' % (filename), 'rb')
        g = open('%s_%s.nc' % (self.params['uuid'], 'dump'), 'wb')
        g.write(f.read())
        f.close()
        g.close()

        return netcdf('%s_%s.nc' % (self.params['uuid'], 'dump'), 'r')


    def close_file(self, data):
        """Closes and deletes temporary file.

        """
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
        if var == None:
            var = self.params['var_list']

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
        
        # Ressets variables
        Var = dict()
        Datasets = dict()
        for item in var:
            Var[item] = ma.zeros(shape)
            try:
                Datasets[self.params['var_tcid'][item][0]][1].append(
                    self.params['var_tcid'][item][2]
                )
            except:
                Datasets[self.params['var_tcid'][item][0]] = [
                    self.params['var_tcid'][item][1],
                    [self.params['var_tcid'][item][2]]
                ]
        
        # Walks through every time index and loads data range from maps.
        for n, T in enumerate(t):
            t2 = time()
            if profile:
                s = '\rLoading data... %s ' % (profiler(shape[0], n + 1, 0, 
                    t1, t2),)
                stdout.write(s)
                stdout.flush()

            # Reads NetCDF file for each dataset
            for Dataset, (Datavar, Datagrid) in Datasets.items():
                params = dict(path=self.params['path'], dataset=Dataset,
                    datavar=Datavar, **self.params['file_list'][N[n]])
                fname = self.create_filename(**params)
                data = self.read_file(fname)
                #
                for Grid in Datagrid:
                    nvar = self.params['var_dict']['{0}_{1}'.format(Dataset,
                        Grid)]
                    if (('lon_i' in self.params.keys()) &
                        ('lat_j' in self.params.keys())):
                        P = data.variables[Grid].data.T[self.params['lat_j'],
                            self.params['lon_i']][JJ, II]
                    else:
                        P = data.variables[Grid].data.T[JJ, II]
                    P[P >= self.variables[item].missing_value] = nan
                    P = ma.masked_where(isnan(P), P)
                    if nonan:
                        P.data[P.mask] = 0
                    #
                    Var[nvar][n, 0, :, :] += P[:, :]
                #
                self.close_file(data)
        
        # If result dictionary contains only one item, return only the value
        # of this item.
        if len(Var.keys()) == 1:
            Var = Var[Var.keys()[0]]
        
        if profile:
            stdout.write('\r\n')
            stdout.flush()
        
        if DEBUG:
            print 't: ', t
            print 'z: ', z
            print 'y:', y
            print 'x:', x
            print 'var: ', Var
            print 'N: ', N
            print 'K: ', K
            print 'J: ', J
            print 'I:', I
            print 'shape: ', shape
        
        if result == 'full':
            return t, z, y, x, Var
        elif result == 'indices':
            return N, K, J, I, Var
        elif result == 'var only':
            return Var
        else:
            raise Warning("Result parameter set imporperly to '%s', "
                "assuming 'var only'." % (result))
            return Var


    def write(self, t, dat):
        """Writes data to NetCDF file.

        PARAMETERS
            t (float) :
                Time.
            dat (dictionary) :
                Data to be saved.
        
        """
        # Some parameters
        attribs = ['amip', 'grib', 'missing_value', 'canonical_units',
            'description', 'standard_name']
        coordinates = dict(k='height', j='latitude', i='longitude')
        # Creates NetCDF file and sets its attributes.
        TM = dates.num2date(t)
        sname = 'tauxy%04d%02d%02d.nc' % (TM.year, TM.month, TM.day)
        f = netcdf('%s/%s' % (self.params['path'], sname), 'w')
        setattr(f, 'name', self.name)
        setattr(f, 'description', self.name)
        for attrib, attrib_value in self.attributes.items():
            setattr(f, attrib, attrib_value)
        # Create a dimension and coordinates of the data.
        f.createDimension('n', 1)
        for dim, coord in coordinates.items():
            if coord == 'n':
                f.createDimension(dim, 1)
                fvar = f.createVariable(coord, 'float', (1, ))
                fvar[:] = t
            else:
                f.createDimension(dim, self.dimensions[dim])
                fvar = f.createVariable(coord, 'float', (dim, ))
                fvar[:] = self.variables[coord].data
            for attrib in attribs:
                try:
                    setattr(fvar, attrib,
                        getattr(self.variables[coord], attrib))
                except:
                    pass
        # Now saves the data
        for var, value in dat.items():
            fvar = f.createVariable(var, 'float', ('n', 'k', 'j', 'i', ))
            for attrib in attribs:
                if getattr(self.variables[var], attrib) != None:
                    setattr(fvar, attrib,
                        getattr(self.variables[var], attrib))
            Value = value.data.copy()
            Value[value.mask] = self.variables[var].missing_value
            fvar[:] = Value[None, None, :, :]
        #
        f.close()
