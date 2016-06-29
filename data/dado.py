# -*- coding: utf-8 -*-
"""Atlantis data framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional data and
associated metadata. Much of this is based uppon Iris library

This module implements a framework for the distributed data archive
DADO running a MySQL database server.

DISCLAIMER
    This software may be used, copied, or redistributed as long as it
    is not sold and this copyright notice is reproduced on each copy
    made. This routine is provided as is without any express or implied
    warranties whatsoever.

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    1 (2013-07-28 12:19 -0300)

"""
from __future__ import division, unicode_literals

__version__ = '$Revision: 1 $'
# $Source$

from time import time
t0 = time() # module start time for profiler

from collections import OrderedDict
from dateutil.parser import parse
from datetime import timedelta
from matplotlib.dates import date2num, epoch2num, num2date
from mysql import connector as mysql
from numpy import asarray, array, fromiter, ma, ndarray, string_, round, float32, float64, isnan
from numpy.core.records import fromarrays
from numpy.ma import MaskedArray, is_masked

import atlantis.data
import atlantis.units

DEBUG = False

###############################################################################
# TODO
#
# createDataset, readDataset
#
# insertData, readData
#
# super class: removeVariable, removeLevel, removeDataset, removeProject,
#              removeStation, removeData

###############################################################################
# PARAMETERS
#

###############################################################################
# CLASSES AND FUNCTIONS
#
class Connection():
    _db = None
    _cursor = None

    def __init__(self, host='localhost', user='root', password='',
        db='dado'):
        self._db = mysql.connect(host=host, user=user, passwd=password, db=db)
        self._cursor = self._db.cursor()
        # This makes sure that conversion of any numpy value class is
        # properly converted according to Stack Overflow user 'mtrbean':
        # http://stackoverflow.com/questions/17053435/
        #   mysql-connector-python-insert-python-variable-to-mysql-table
        self._db.set_converter_class(NumpyMySQLConverter)
        # Important!! changes time zone to UTC
        self.execute('SET time_zone=\'+00:00\'')


    def close(self):
        self._db.close()


    def escape_none(self, var):
        """Escape null variables."""
        if var == None:
            return "NULL"
        return var


    def _execute(self, sql, data=(), fetch=True, commit=False, debug=False,
        **kwargs):
        if isinstance(data, list):
            self._cursor.executemany(sql, data, **kwargs)
        else:
            self._cursor.execute(sql, data, **kwargs)
        if debug:
            print(self._cursor._executed)
        if commit:
            try:
                self._db.commit()
            except:
                pass
        if fetch:
            return self._cursor.fetchall()
        else:
            return None


    def execute(self, sql, data=(), commit=True, debug=False):
        self._execute(sql, data, commit=commit, fetch=False, debug=debug)


    def fetch(self, sql, data=(), commit=False, debug=False):
        return self._execute(sql, data, commit=commit, fetch=True, debug=debug)


    def datavars_read(self, name=None, where=None, mode='list', debug=False):
        """Reads list of data variables.

        Parameters
        ----------
        name : string or sequence, optional
            String or sequence of strings with the names of the
            variables to retrieve.
        where : string or sequence, optional
            List of SQL-like selection arguments. If given a string, it
            will become the SQL argument. If given a list of arguments,
            they will be joined using 'OR'.
        mode : string, optional
            If `list` (default), returns a list of variables. If
            `dict`, returns a dictionary with variables with variable
            names as keys.

        debug : bool, optional

        Returns
        -------
        result : sequence
            A single variable, a list or dictionary of variables.

        Examples
        --------
        var_list = datavars_read(where=['var_id=10010'])
        var_list = datavars_read(name='eastward_wind'])

        """
        _sql = 'SELECT * FROM datavars '
        if isinstance(where, list) | isinstance(where, tuple):
            where = ' OR '.join(where)
        elif where is None:
            where = ''
        if isinstance(name, list) | isinstance(name, tuple):
            where += ' OR '.join(['name="{}"'.format(_name) for _name in name])
        elif isinstance(name, basestring):
            if where != '':
                where += ' OR '
            where += 'name="{}"'.format(name)
        if isinstance(where, basestring):
            _sql += 'WHERE ({});'.format(where)

        # Remember that the sequence of columns is as follows:
        # var_id, name, canonical_units, units, scale_factor, add_offset,
        # long_name, symbol, description, dictionary
        _result = list()
        for _item in self.fetch(_sql, debug=debug):
            _var = atlantis.data.Variable(
                id = _item[0],
                standard_name = _item[1],
                canonical_units = _item[2],
                units = _item[3],
                scale_factor = _item[4],
                add_offset = _item[5],
                long_name = _item[6],
                symbol = _item[7],
                description = _item[8]
            )
            _result.append(_var)

        if mode == 'list':
            if len(_result) == 0:
                return None
            elif len(_result) == 1:
                return _result[0]
            else:
                return _result
        elif mode == 'dict':
            return OrderedDict([(_item.standard_name, _item)
                for _item in _result])
        else:
            raise ValueError('Invalid mode `{}`.'.format(mode))


    def datavars_nextId(self, var_name=None):
        """Returns next available data variable ID."""
        _sql = 'SELECT max(var_id) as max_var_id FROM dado.datavars;'
        try:
            _result = self.fetch(_sql)[0][0] + 1
        except:
            _result = 10001
        return _result


    def datavars_insert(self, variables, update=True):
        """
        Inserts and updates data variables to database.

        Parameters
        ----------
        variables : Variable or list
            Creates a single or a list of variables.
        update : boolean, optional
            If true (default), updates existing variables.

        Returns
        -------
        Nothing

        """
        if not isinstance(variables, list):
            variables = [variables]

        for _item in variables:
            # Checks if variable ID has been set. If not, creates one.
            if _item.id == None:
                _item.id = self.datavars_nextId()
            _sql = ('INSERT INTO datavars '
                '(var_id, name, canonical_units, units, scale_factor, '
                ' add_offset, long_name, symbol, description) '
                'VALUES(%(id)s, %(standard_name)s, %(canonical_units)s, '
                '%(units)s, %(scale_factor)s, %(add_offset)s, %(long_name)s, '
                '%(symbol)s, %(description)s)'
            )
            _data = _item.__dict__
            if update:
                _sql += (' ON DUPLICATE KEY UPDATE '
                    #'var_id=VALUES(var_id), '
                    'name=VALUES(name), '
                    'canonical_units=VALUES(canonical_units), '
                    'units=VALUES(units), '
                    'scale_factor=VALUES(scale_factor), '
                    'add_offset=VALUES(add_offset), '
                    'long_name=VALUES(long_name), '
                    'symbol=VALUES(symbol), '
                    'description=VALUES(description)'
                )
            _sql += ';'
            self.execute(_sql, _data, commit=True)


    def datalevels_insert(self, levels, update=True):
        """
        Inserts new levels to database.

        PARAMETERS
            levels (Level or list):
                Creates a single of a list of levels.
            update (boolean, optional) :
                If true (default), updates existing levels.

        RESULTS
            Nothing

        """
        if not isinstance(levels, list):
            levels = [levels]

        for _item in levels:
            _sql = ('INSERT INTO datalevels '
                '(level_id, name, temporal, horizontal_x, horizontal_y, '
                'vertical, description) '
                'VALUES(%(id)s, %(name)s, %(temporal)s, %(horizontal_x)s, '
                '%(horizontal_y)s, %(vertical)s, %(description)s)'
            )
            if update:
                _sql += (' ON DUPLICATE KEY UPDATE '
                    'name=VALUES(name), '
                    'temporal=VALUES(temporal), '
                    'horizontal_x=VALUES(horizontal_x), '
                    'horizontal_y=VALUES(horizontal_y), '
                    'vertical=VALUES(vertical), '
                    'description=VALUES(description)'
                )
            _sql += ';'
            _data = _item.attributes
            self.execute(_sql, _data)


    def datalevels_read(self, level=None, where=None, debug=False):
        """Reads list of data levels.

        PARAMETERS
            where (string or list, optional) :
                List of SQL-like selection arguments. If given a string,
                it will become the SQL argument. If given a list of
                arguments, they will be joined using 'OR'.

        RETURNS
            result (list) :
                A single level or a list of levels.

        EXAMPLES:
            level_list = datalevels_read(where=['level_id=0'])
            level_list = datalevels_read(where=['name='raw'])

        """
        _sql = 'SELECT * FROM datalevels '
        if isinstance(where, list):
            where = ' OR '.join(where)
        if isinstance(where, basestring):
            _sql += 'WHERE ({});'.format(where)
        elif isinstance(level, basestring):
            _sql += 'WHERE (name=\'{}\');'.format(level)

        # Remember that the sequence of columns is as follows:
        # (level_id, name, temporal, horizontal_x, horizontal_y, vertical,
        # description, dictionary)
        _result = list()
        for _item in self.fetch(_sql, debug=debug):
            _var = Level(
                id = _item[0],
                name = _item[1],
                temporal = _item[2],
                horizontal_x = _item[3],
                horizontal_y = _item[4],
                vertical = _item[5],
                description = _item[6]
            )
            _result.append(_var)

        if len(_result) == 0:
            return None
        elif len(_result) == 1:
            return _result[0]
        else:
            return _result


    def projects_insert(self, projects, update=True):
        """
        Inserts new projects to database.

        PARAMETERS
            projects (Project or list):
                Creates a single of a list of projects.
            update (boolean, optional) :
                If true (default), updates existing projects.

        RESULTS
            Nothing

        """
        if not isinstance(projects, list):
            projects = [projects]

        for _item in projects:
            _sql = ('INSERT INTO projects '
                '(project_id, parent_id, name, description) '
                'VALUES(%(id)s, %(parent)s, %(name)s, %(description)s)'
            )
            if update:
                _sql += (' ON DUPLICATE KEY UPDATE '
                    'parent_id=VALUES(parent_id), '
                    'name=VALUES(name), '
                    'description=VALUES(description)'
                )
            _sql += ';'
            _data = _item.attributes
            self.execute(_sql, _data)


    def projects_read(self, project=None, where=None, debug=False):
        """
        Reads list of projects.

        PARAMETERS
            where (string or list, optional) :
                List of SQL-like selection arguments. If given a string,
                it will become the SQL argument. If given a list of
                arguments, they will be joined using 'OR'.

        RETURNS
            result (list) :
                A single project or a list of projects.

        EXAMPLES:
            project_list = projects_read(where=['project_id=0'])
            project_list = projects_read(where=['name='dummy',
                'project_id=null'])

        """
        _result = list()
        #
        if isinstance(project, int) & (where == None):
            where = 'project_id={}'.format(project)
        #
        if where != None:
            _sql = 'SELECT * FROM projects '
            if isinstance(where, list) | isinstance(where, tuple):
                where = ' OR '.join(where)
            if isinstance(where, basestring):
                _sql += 'WHERE ({});'.format(where)
            # Remember that the sequence of columns is as follows:
            # (project_id, parent_id, name, description, dictionary)
            for _item in self.fetch(_sql, debug=debug):
                _var = Project(
                    id=_item[0],
                    parent=_item[1],
                    name=_item[2],
                    description=_item[3]
                )
                _result.append(_var)
        else:
            if project == None:
                project = []
            elif isinstance(project, basestring):
                project = [project]
            #
            for item in project:
                url = item.split('/')
                parent_id = None
                # Walk through every parent project
                _var = None
                for sub_item in url:
                    if parent_id in ['NULL', None]:
                        _sql = ('SELECT * FROM projects WHERE (name="{}" AND '
                            'parent_id IS NULL)'.format(sub_item))
                    else:
                        _sql = ('SELECT * FROM projects WHERE (name="{}" AND '
                            'parent_id={})'.format(sub_item, parent_id))
                    # Remember that the sequence of columns is as follows:
                    # (project_id, parent_id, name, description, dictionary)
                    fetch = self.fetch(_sql, debug=debug)[0]
                    parent_id = fetch[0] # Current project ID is new parent ID!
                    _var = Project(
                        id = fetch[0],
                        parent = _var,
                        name = fetch[2],
                        description = fetch[3]
                    )
                #
                _result.append(_var)

        if len(_result) == 0:
            return None
        elif len(_result) == 1:
            return _result[0]
        else:
            return _result


    def stations_insert(self, stations, update=True):
        """
        Inserts new stations to database.

        PARAMETERS
            stations (Station or list):
                Creates a single of a list of stations.
            update (boolean, optional) :
                If true (default), updates existing stations.

        RESULTS
            Nothing

        """
        if not isinstance(stations, list):
            stations = [stations]

        for _item in stations:
            _sql = ('INSERT INTO stations '
                '(project_id, station_id, name, x, y, z, description) '
                'VALUES(%(project)s, %(id)s, %(name)s, %(x)s, %(y)s, '
                '%(z)s, %(description)s)'
            )
            if update:
                _sql += (' ON DUPLICATE KEY UPDATE '
                    'name=VALUES(name), '
                    'x=VALUES(x), '
                    'y=VALUES(y), '
                    'z=VALUES(z), '
                    'description=VALUES(description)'
                )
            _sql += ';'
            # Rounds longitude and latitude coordinates to 8 decimal places.
            try:
                _item.x = round(_item.x, 8)
            except:
                _item.x = None
            try:
                _item.y = round(_item.y, 8)
            except:
                _item.y = None
            _data = dict(
                project = _item.project.id,
                id = _item.id,
                name = _item.name,
                x = _item.x,
                y = _item.y,
                z = _item.z,
                description = _item.description
            )
            self.execute(_sql, _data)


    def stations_read(self, station=None, project=None, where=None,
        create=None, result='list', debug=False):
        """Reads list of stations.

        PARAMETERS
            station (string) :
                Name of the station.
            project (Project, optional) :
                Project of the station.
            where (string or list, optional) :
                List of SQL-like selection arguments. If given a string,
                it will become the SQL argument. If given a list of
                arguments, they will be joined using 'OR'.
            create (Station, optional) :
                If given and query results are empty, then creates
                station.
            result (string, optional) :
                If set to `list`, returns a list of stations, if set to
                `dict`, returns a dictionary of stations with
                `project.id`.`station.id` as keys.

        RETURNS
            result (list) :
                A single station or a list of stations.

        EXAMPLES:
            stations_list = stations_read(where=['stations_id=0'])
            stations_list = stations_read(where=['name='free',
                'project_id=null'])

        """
        # TODO: station and project
        _where = ''
        if project != None:
            if isinstance(project, Project):
                _where += 'project_id={}'.format(project.id)
            elif isinstance(project, int):
                _where += 'project_id={}'.format(project)
            elif isinstance(project, basestring):
                project = self.projects_read(project)
                _where += 'project_id={}'.format(project.id)
            else:
                raise ValueError('Invalid project `{}`.'.format(project))
        if station != None:
            if _where != '':
                _where += ' AND '
            _where += 'name=\'{}\''.format(station)
        #
        if isinstance(where, list) == list:
            where = ' OR '.join(where)
        if isinstance(where, basestring):
            if _where != '':
                _where += '({}) AND ({})'.format(_where, where)
            else:
                _where += '{}'.format(where)
        #
        if _where == '':
            return None
        else:
            _sql = ('SELECT * FROM stations WHERE {} ORDER BY project_id, '
                'station_id;').format(_where)

        # Remember that the sequence of columns is as follows:
        # (project_id, station_id, name, x, y, z, description, dictionary)
        _result = list()
        _ids = list()
        for _item in self.fetch(_sql, debug=debug):
            _var = Station(
                project = self.projects_read(project=_item[0], debug=debug),
                id = _item[1],
                name = _item[2],
                x = _item[3],
                y = _item[4],
                z = _item[5],
                description = _item[6]
            )
            _result.append(_var)
            _ids.append('{}.{}'.format(_var.project.id, _var.id))

        if len(_result) == 0:
            if create != None:
                # Creates the station and returns it.
                create.id = self.stations_nextId(project=create.project)
                self.stations_insert(create)
                return create
            else:
                return None
        elif (len(_result) == 1) & (result == 'list'):
            return _result[0]
        else:
            if result == 'dict':
                _result = OrderedDict(zip(_ids, _result))
            return _result


    def stations_nextId(self, project, debug=False):
        """Returns next available station ID for given project."""
        _sql = ('SELECT max(station_id) as max_station_id FROM '
            'dado.stations WHERE project_id={};').format(project.id)
        try:
            _result = self.fetch(_sql, debug=debug)[0][0] + 1
        except:
            _result = 1001
        return _result


    def instruments_read(self, where=None, result='list', debug=False):
        """Reads list of instruments.

        PARAMETERS
            where (string or list, optional) :
                List of SQL-like selection arguments. If given a string,
                it will become the SQL argument. If given a list of
                arguments, they will be joined using 'OR'.
            result (string, optional) :
                If `list`, returns a list of instruments. If `dict`,
                returns a dictionary of instruments with instrument ID
                as keys.

        RETURNS
            result (list or dictionary) :
                A single instrument or a list of instruments.

        EXAMPLES:
            instruments_list = instruments_read()
            instruments_list = instruments_read(where='product_name='C3')

        """
        _sql = 'SELECT * FROM instruments '
        if isinstance(where, list):
            where = ' OR '.join(where)
        if isinstance(where, basestring):
            _sql += 'WHERE ({});'.format(where)

        # Remember that the sequence of columns is as follows:
        # (project_id, instrument_id, manufacturer, product_name,
        # serial_number, institution, description)
        if result == 'list':
            _result = list()
        elif result == 'dict':
            _result = dict()
        else:
            raise ValueError('Invalid return parameter `{}`.'.format(result))
        for _item in self.fetch(_sql, debug=debug):
            _var = Instrument(
                project = self.projects_read(project=_item[0]),
                id = _item[1],
                manufacturer = _item[2],
                name = _item[3],
                serial = _item[4],
                institution = _item[5],
                description = _item[6]
            )
            if result == 'list':
                _result.append(_var)
            elif result == 'dict':
                _result[_var.id] = _var

        if len(_result) == 0:
            return None
        elif len(_result) == 1:
            return _result[0]
        else:
            return _result


    def dataset_insert(self, dataset, update=True):
        """
        Inserts new dataset to database.

        PARAMETERS
            dataset (Level or list):
                The dataset object.
            update (boolean, optional) :
                If true (default), updates existing dataset.

        RESULTS
            Nothing

        """
        # First step: add items to `datasets` table:
        _sql = (
            'INSERT INTO datasets '
            '(dataset_id, name, description) '
            'VALUES(%(id)s, %(name)s, %(description))'
        )
        if update:
            _sql += (
                'ON DUPLICATE KEY UPDATE '
                'name=VALUES(name), '
                'description=VALUES(description)'
            )

        _sql += ';'
        _data = dataset.attributes
        self.execute(_sql, _data)

        # Second step: add items to `datasetvars` table:
        raise Warning('This is not implemented yet!.')
        for _item in dataset.var_list:
            _sql = ('INSERT INTO datasevars '
                '(project_id, level_id, var_id) '
                'VALUES(%(project)s, %(level)s, %(var)s)'
            )
            if update:
                _sql += (' ON DUPLICATE KEY UPDATE '
                    'name=VALUES(name), '
                    'temporal=VALUES(temporal), '
                    'horizontal_x=VALUES(horizontal_x), '
                    'horizontal_y=VALUES(horizontal_y), '
                    'vertical=VALUES(vertical), '
                    'description=VALUES(description)'
                )
            _sql += ';'
            _data = _item.attributes
            self.execute(_sql, _data)


class Level(object):
    """
    Data level class.

    ATTRIBUTES
        id           -- Identification of the level in database.
        name         -- Name of the level.
        temporal     -- Temporal resolution (i.e. seconds, hours, days)
        horizontal_x -- Zonal horizontal resolution (i.e. meters,
                        degrees)
        horizontal_y -- Meridional horizontal resolution (i.e. meters,
                        degrees)
        vertical     -- Vertical resolution in meters.
        description  -- Long description of level.

    According to the use, spatial resolution should be either given in
    meters or degrees. To ensure consistency with matplotlib's date
    format, days should be used as time scale.

    """
    def __init__(self, **kwargs):
        self.attributes = dict()
        #
        for key, item in kwargs.items():
            try:
                setattr(self, key, item)
            except:
                print 'Warning: Invalid attribute {0}'.format(key)
                pass
        #
        return


    def _set_attribute(self, attrib, val):
        self.attributes[attrib] = val


    def _get_attribute(self, attrib):
        if attrib in self.attributes.keys():
            return self.attributes[attrib]
        else:
            return None


    @property
    def id(self):
        """Id of the level."""
        return self._get_attribute('id')
    @id.setter
    def id(self, val):
        self._set_attribute('id', val)

    @property
    def name(self):
        """Name of the level."""
        return self._get_attribute('name')
    @name.setter
    def name(self, val):
        self._set_attribute('name', val)

    @property
    def temporal(self):
        """Temporal resolution."""
        return self._get_attribute('temporal')
    @temporal.setter
    def temporal(self, val):
        self._set_attribute('temporal', val)

    @property
    def horizontal_x(self):
        """Zonal horizontal resolution."""
        return self._get_attribute('horizontal_x')
    @horizontal_x.setter
    def horizontal_x(self, val):
        self._set_attribute('horizontal_x', val)

    @property
    def horizontal_y(self):
        """Meridional horizontal resolution."""
        return self._get_attribute('horizontal_y')
    @temporal.setter
    def horizontal_y(self, val):
        self._set_attribute('horizontal_y', val)

    @property
    def vertical(self):
        """Vertical resolution."""
        return self._get_attribute('vertical')
    @vertical.setter
    def vertical(self, val):
        self._set_attribute('vertical', val)

    @property
    def description(self):
        """Level description."""
        return self._get_attribute('description')
    @description.setter
    def description(self, val):
        self._set_attribute('description', val)


class Project(object):
    """Data project class."""
    def __init__(self, **kwargs):
        self.attributes = dict()
        #
        for key, item in kwargs.items():
            try:
                setattr(self, key, item)
            except:
                print 'Warning: Invalid attribute {0}'.format(key)
                pass
        #
        return


    def _set_attribute(self, attrib, val):
        self.attributes[attrib] = val


    def _get_attribute(self, attrib):
        if attrib in self.attributes.keys():
            return self.attributes[attrib]
        else:
            return None


    @property
    def id(self):
        """The id of the project."""
        return self._get_attribute('id')
    @id.setter
    def id(self, val):
        self._set_attribute('id', val)

    @property
    def parent(self):
        """The parent of the project."""
        return self._get_attribute('parent')
    @parent.setter
    def parent(self, val):
        self._set_attribute('parent', val)

    @property
    def name(self):
        """The name of the project."""
        return self._get_attribute('name')
    @name.setter
    def name(self, val):
        self._set_attribute('name', val)

    @property
    def description(self):
        """The description of the project."""
        return self._get_attribute('description')
    @description.setter
    def description(self, val):
        self._set_attribute('description', val)


class Station(object):
    """Station class."""
    def __init__(self, **kwargs):
        self.attributes = dict()
        #
        for key, item in kwargs.items():
            try:
                setattr(self, key, item)
            except:
                print 'Warning: Invalid attribute {0}'.format(key)
                pass
        #
        return


    def _set_attribute(self, attrib, val):
        self.attributes[attrib] = val


    def _get_attribute(self, attrib):
        if attrib in self.attributes.keys():
            return self.attributes[attrib]
        else:
            return None



    @property
    def project(self):
        """The project of the station."""
        return self._get_attribute('project')
    @project.setter
    def project(self, val):
        self._set_attribute('project', val)

    @property
    def id(self):
        """The id of the project."""
        return self._get_attribute('id')
    @id.setter
    def id(self, val):
        self._set_attribute('id', val)

    @property
    def name(self):
        """The name of the project."""
        return self._get_attribute('name')
    @name.setter
    def name(self, val):
        self._set_attribute('name', val)

    @property
    def x(self):
        """Zonal coordinate."""
        return self._get_attribute('x')
    @x.setter
    def x(self, val):
        self._set_attribute('x', val)

    @property
    def y(self):
        """Meridional coordinate."""
        return self._get_attribute('y')
    @y.setter
    def y(self, val):
        self._set_attribute('y', val)

    @property
    def z(self):
        """Depth."""
        return self._get_attribute('z')
    @z.setter
    def z(self, val):
        self._set_attribute('z', val)

    @property
    def description(self):
        """The description of the project."""
        return self._get_attribute('description')
    @description.setter
    def description(self, val):
        self._set_attribute('description', val)


class Instrument(object):
    """Instrument class."""
    def __init__(self, **kwargs):
        self.attributes = dict()
        #
        for key, item in kwargs.items():
            try:
                setattr(self, key, item)
            except:
                print 'Warning: Invalid attribute {0}'.format(key)
                pass
        #
        return


    def _set_attribute(self, attrib, val):
        self.attributes[attrib] = val


    def _get_attribute(self, attrib):
        if attrib in self.attributes.keys():
            return self.attributes[attrib]
        else:
            return None


    @property
    def project(self):
        """The project of the instrument."""
        return self._get_attribute('project')
    @project.setter
    def project(self, val):
        self._set_attribute('project', val)

    @property
    def id(self):
        """The ID of the instrument."""
        return self._get_attribute('id')
    @id.setter
    def id(self, val):
        self._set_attribute('id', val)

    @property
    def manufacturer(self):
        """The manufacturer of the instrument."""
        return self._get_attribute('manufacturer')
    @manufacturer.setter
    def manufacturer(self, val):
        self._set_attribute('manufacturer', val)

    @property
    def name(self):
        """The name of the instrument."""
        return self._get_attribute('name')
    @name.setter
    def name(self, val):
        self._set_attribute('name', val)

    @property
    def serial(self):
        """The serial number of the instrument."""
        return self._get_attribute('serial')
    @serial.setter
    def serial(self, val):
        self._set_attribute('serial', val)

    @property
    def institution(self):
        """The institution, owner of the instrument."""
        return self._get_attribute('institution')
    @institution.setter
    def institution(self, val):
        self._set_attribute('institution', val)

    @property
    def description(self):
        """The description of the instrument."""
        return self._get_attribute('description')
    @description.setter
    def description(self, val):
        self._set_attribute('description', val)


class Dataset(object):
    """Dataset class."""
    def __init__(self, **kwargs):
        self.attributes = dict()
        #
        for key, item in kwargs.items():
            try:
                setattr(self, key, item)
            except:
                print 'Warning: Invalid attribute {0}'.format(key)
                pass
        #
        return


    def add_variable(self, var):
        """Adds a variable to the dataset."""
        var_list = self.var_list
        if not isinstance(var_list, dict):
            var_list = dict()
        var_list[var.standard_name] = var
        self.var_list = var_list


    def get_variable(self, name=None):
        """Returns the variable dictionary or a single variable."""
        if name == None:
            return self.var_list
        elif name in self.var_list.keys():
            return self.var_list[name]
        else:
            return None


    def _set_attribute(self, attrib, val):
        self.attributes[attrib] = val


    def _get_attribute(self, attrib):
        if attrib in self.attributes.keys():
            return self.attributes[attrib]
        else:
            return None


    @property
    def id(self):
        """The dataset ID."""
        return self._get_attribute('id')
    @id.setter
    def id(self, val):
        self._set_attribute('id', val)

    @property
    def name(self):
        """A short name for the dataset."""
        return self._get_attribute('name')
    @name.setter
    def name(self, val):
        self._set_attribute('name', val)

    @property
    def description(self):
        """The description of the dataset."""
        return self._get_attribute('description')
    @description.setter
    def description(self, val):
        self._set_attribute('description', val)

    @property
    def var_list(self):
        """List of variables in dataset."""
        return self._get_attribute('var_list')
    @var_list.setter
    def var_list(self, val):
        self._set_attribute('var_list', val)


class NumpyMySQLConverter(mysql.conversion.MySQLConverter):
    """A mysql.connector converter that handles Numpy types."""

    #def __init__(self, charset=None, use_unicode=True):
    #    mysql.conversion.MySQLConverter.__init__(self)
    #    self.python_types[float32] = self._float64_to_mysql
    #    self.python_types[float64] = self._float64_to_mysql
    #    self.python_types[string_] = self._string_to_mysql

    def _float32_to_mysql(self, value):
        return float(value)

    def _float64_to_mysql(self, value):
        return float(value)

    def _int32_to_mysql(self, value):
        return int(value)

    def _int64_to_mysql(self, value):
        return int(value)

    def _string__to_mysql(self, value):
        return str(value)

   # def _DATETIME_to_python(self, value, desc=None):
   #     """Returns DATETIME column type as matplotlib date number."""
   #     return value

    def _DECIMAL_to_python(self, value, desc=None):
        """Returns value as a floating point number."""
        return float(value)


class Sequence(atlantis.data.Sequence):
    """Common sequencial data for DADO data server.

    """
    _conn = None
    def __init__(self, connection):
        """
        Initializes the sequence class for reading DADO data server
        contents.

        PARAMETERS
            connection (Connection) :
                Connection to the database.

        """
        self._conn = connection


    def write(self, t, z, y, x, val, station, level, var, update=True,
        debug=False, **kwargs):
        """
        Writes data sequence. Input arrays have to be
        one-dimensional. It assumes that all input data is from the
        same station and same level.

        Note that height is positive above sea level and negative below.
        Longitude and latitude coordinates are rounded to the eigth decimal
        place.

        Parameters
        ----------
        t, z, y, x (array like) :
            Time, height, latitude, longitude of the data.
        val (array like) :
            Data value.
        station (Station) :
            Project station object.
        level (Level) :
            Data level.
        var (Variable) :
            Variable of the data.
        update (boolean, optional) :
            If true (default), updates existing stations.
        ms (integer, optional) :
            Milliseconds added add to time.
        minimum, maximum, stdev (float, optional) :
        quality (array like, optional) :
        source (int, optional) :
        instrument (int, optional) :

        Returns
        -------
        Nothing.

        """
        # TODO: There is some kind of inconsistency in saving the data,
        # minimum, maximum and standard deviation values. They should be
        # passed over a Variable instance and not a value.

        # Defines default optional arguments and merges them with delivered
        # optinal arguments
        _default_args = dict(
            ms = 0,
            minimum = None,
            maximum = None,
            stdev = None,
            quality = None,
            source = None,
            instrument = None
        )
        kwargs = dict(_default_args.items() + kwargs.items())

        # This is to make sure that all coordinates arrays are of the same
        # type and have the same dimensions. If single data values are used,
        # it assumes the same values are repeated over the whole dataset.
        # Latitude and longitude coordinates are rounded to eigth decimal
        # digits.
        t = self._asarray(t, dtype='time')
        N = t.size
        z = self._asarray(z, size=N)
        y = self._asarray(y, size=N)
        try:
            y = round(y, 8)
        except:
            pass
        x = self._asarray(x, size=N)
        try:
            x = round(x, 8)
        except:
            pass
        try:
            z = round(z, 4)
        except:
            pass
        if (z.size != N) | (y.size != N) | (x.size != N):
            raise ValueError('Size of one of the arrays does not match.')

        # Makes sure that level is of Level class
        if isinstance(level, basestring):
            level = self._conn.datalevels_read(level, debug=debug)

        # In case multiple variables are inserted at the same time, converts
        # all input to the same format: a dictionary containing the variable
        # standard name as key and a structured array for data values. First
        # treats the variables and then the values.
        if isinstance(var, atlantis.data.Variable):
            var = OrderedDict([(var.standard_name, var)])
        elif isinstance(var, list) | isinstance(var, tuple):
            var = OrderedDict([(_item.standard_name, _item) for _item in var])
        #
        nvars = len(var)
        val = self._asarray(val)
        if nvars == 1:
            val = val[None, :]
        val = fromarrays(val,
            dtype={'names': var.keys(), 'formats':['f4']*nvars})

        # Makes sure data, minimum, maximum and standard deviation are saved
        # in canonical units.
        for key in var.keys():
            if var[key].units == None:
                _from = var[key].canonical_units
            else:
                _from = var[key].units
            _to = var[key].canonical_units
            val[key] = atlantis.units.fromAtoB(val[key], _from, _to)
        for key in ['minimum', 'maximum', 'stdev']:
            if kwargs[key] != None:
                kwargs[key] = atlantis.units.fromAtoB(kwargs[key], _from, _to,
                vtype=key)

        # Walks through each coordinate and writes data from each variable
        # to database
        _sql = ('INSERT INTO datavalues '
            '(project_id, station_id, level_id, var_id, x, y, z, t, '
            'ms, minimum, value, maximum, stdev, quality, source, '
            'instrument, modified_on) '
            #
            'VALUES(%(project)s, %(station)s, %(level)s, %(var)s, '
            '%(x)s, %(y)s, %(z)s, %(t)s, %(ms)s, %(min_)s, %(val)s, '
            '%(max_)s, %(stdev)s, %(quality)s, %(source)s, '
            '%(instrument)s, NOW())'
        )
        _data = []
        for i, (tt, zz, yy, xx) in enumerate(zip(t, z, y, x)):
            for vv in var.values():
                _val = self._read_asarray(val[vv.standard_name], i)
                if _val == None:
                    continue
                _data.append(dict(
                    project = station.project.id,
                    station = station.id,
                    level = level.id,
                    var = vv.id,
                    x = xx,
                    y = yy,
                    z = zz,
                    t = tt,
                    ms = self._read_asarray(kwargs['ms'], i),
                    min_ = self._read_asarray(kwargs['minimum'], i),
                    val = _val,
                    max_ = self._read_asarray(kwargs['maximum'], i),
                    stdev = self._read_asarray(kwargs['stdev'], i),
                    quality = self._read_asarray(kwargs['quality'], i),
                    source = self._read_asarray(kwargs['source'], i),
                    instrument = self._read_asarray(kwargs['instrument'], i)
                ))
        if update:
            _sql += (' ON DUPLICATE KEY UPDATE '
                'minimum=VALUES(minimum), '
                'value=VALUES(value), '
                'maximum=VALUES(maximum), '
                'stdev=VALUES(stdev), '
                'quality=VALUES(quality), '
                'source=VALUES(source), '
                'instrument=VALUES(instrument), '
                'modified_on=NOW()'
        )
        _sql += ';'
        self._conn.execute(_sql, _data, debug=debug)
        # That's it!


    def update(self, t, z, y, x, station, level, var, new_val, debug=False,
        **kwargs):
        """
        Updates database entry. Input arrays have to be one-dimensional.
        It assumes that all input data is from the same station and same
        level.

        Note that height is positive above sea level and negative below.
        Longitude and latitude coordinates are rounded to the eigth
        decimal place.

        Parameters
        ----------
        t, z, y, x : array like
            Time, height, latitude, longitude of the data.
        station : Station
            Project station object.
        level : Level
            Data level.
        var : Variable
            Variable of the data.
        new_val : array like, dictionary of arrays.
            New updated data value. If an array is given, then updates
            the value for each entry. If a dictionary of arrays is
            given, the fields to be updated are the keys to the
            dictionary.

        Returns
        -------
        Nothing.

        """
        # This is to make sure that all coordinates arrays are of the same
        # type and have the same dimensions. If single data values are used,
        # it assumes the same values are repeated over the whole dataset.
        # Latitude and longitude coordinates are rounded to eigth decimal
        # digits.
        t = self._asarray(t, dtype='time')
        N = t.size
        z = self._asarray(z, size=N)
        y = self._asarray(y, size=N)
        try:
            y = round(y, 8)
        except:
            pass
        x = self._asarray(x, size=N)
        try:
            x = round(x, 8)
        except:
            pass
        if (z.size != N) | (y.size != N) | (x.size != N):
            raise ValueError('Size of one of the arrays does not match.')

        # Defines default optional arguments and merges them with delivered
        # optinal arguments
        _default_args = dict(
            ms = 0
        )
        kwargs = dict(_default_args.items() + kwargs.items())

        # Checks whether `new_val` is a dictionary.
        if not isinstance(new_val, dict):
            new_val = dict(value=new_val)

        # Builds default SQL strings.
        _update = 'UPDATE dado.datavalues SET '
        _where = (' WHERE project_id={} AND station_id={} AND level_id={} AND '
            'var_id={} AND x={} AND y={} AND z={} AND t="{}" AND ')

        # Important data keys
        update_keys = new_val.keys()
        where_keys = kwargs.keys()

        # Walks through each entry in new value dictionary and checks the size
        for i in xrange(N):
            _set = ', '.join(['{}={}'.format(key,
                self._read_asarray(new_val[key], i)) for key in update_keys])
            _erehw = ' AND '.join(['{}={}'.format(key,
                self._read_asarray(kwargs[key], i)) for key in where_keys])
            _sql = (
                _update +
                _set +
                _where.format(station.project.id, station.id, level.id, var.id,
                    x[i], y[i], z[i], t[i].strftime('%Y-%m-%d %H:%M:%S')) +
                _erehw +
                ';'
            )
            _sql = _sql.replace('=None', ' IS NULL')
            self._conn.execute(_sql, debug=debug)
        #
        # That's it!


    def read(self, var, t=None, z=None, y=None, x=None,
        radius_t=0., radius_v=0., radius_h=0.,
        tlim=None, zlim=None, ylim=None, xlim=None,
        level=None, project=None, station=None, where=None,
        column='value', skip=[], units='units', sort=True,
        groupby='t, ms, z, y, x', quality=1, masked=False, profile=True,
        debug=False):
        """
        Reads sequence in dataset.

        Parameters
        ----------
        var : string, array like, Variable
            Variables to be read from dataset. If a single
            string or a list of strings are given, uses them as
            variable names. If a single or list of class
            Variable is given, uses appropriate variable IDs.
        t, z, y, z : array like, optional
            List of temporal, vertical, meridional and zonal
            point coordinates of interest.
        radius_t, radius_v, radius_h : float, optional
            Temporal, vertical and horizontal search radii in
            days, meters and degrees, respectively.
        tlim, zlim, ylim, xlim : sequence , optional
            The temporal, vertical, meridional and zonal limits
            (minimum, maximum) for which data will be read.
        level : string, array like, Level, optional
            Single level of list of levels to read data from. If
            a single string of list of strings are given, uses
            them as level names. If a single or list of class
            Level is given, uses appropriate level IDs.
        project : string, array like, Project, optional
            Single project of list of projects to read data from.
            If a single string of list of strings are given,
            uses them as project names. If a single or list
            of class Project is given, uses appropriate
            variable IDs.
        station : string, array like, Station, optional
            Single station of list of stations to read data from.
            If a single string of list of strings are given,
            uses them as station names. If a single or list
            of class Station is given, uses appropriate
            project and station IDs.
        where : string, optional
        column : string, optional
            Indicates from which database column data will be
            read: minimum, value (default), maximum, stdev,
            quality, source, instrument.
        skip : array like, optional
            Skips data column.
        units : string, list, dictionary, optional
            Gives the desired units for the variables to be read.
            Accepts either `canonical_units`, `units`, a list or
            dictionary of units. If a list is geven, items should be
            in the same order as `var`. If it is a dictionary, keys
            should be the same as in `var` as keys. If not given,
            converts canonical units to units.
        sort : boolean, string optional
            If true, sorts the data record in order of ascendant
            time, latitude, longitude and height (depth). Also allows
            to give a string sequence of columns separated by commas.
        quality : int, optional
            Gives the quality flag to read.
        masked : boolean, optional
            If true, returns masked array and masks values with NaN.
        profile : boolean, optional
            Sets whether the status is send to screen.
        debug : boolean, optional
            If true, shows SQL commands on screen.

        Returns
        -------
        dat : structured array
            Record time-series of 'time', 'latitude', 'longitude',
            selected variable and 'mission'.
        variables : OrderedDict
            Variable definitions.

        """
        # Checks project and station parameters.
        if isinstance(project, basestring):
            project = self._conn.projects_read(project=project, debug=debug)
        if isinstance(station, basestring):
            station = self._conn.stations_read(project=project,
                station=station, debug=debug)

        if column not in ['minimum', 'value', 'maximum', 'stdev', 'quality',
            'source', 'instrument']:
            raise ValueError('Invalid column `{}`'.format(column))
        # Checks parameter `var` and, if appropriate, converts it to class
        # Variable. Afterwards generetas `where` and pivot SQL clauses.
        _where = ''
        _pivots = ''
        _var_return = OrderedDict()
        _units = dict()
        # Defines the type of data columns. Note that instrument and quality
        # are set to `float64` since they can have null (`None`) values which
        # will be converted no `nan`.
        columns = ['project_id', 'station_id', 'x', 'y', 'z', 't', 'ms',
            'instrument', 'quality']
        _var_list, _dtype, _columns = [], [], []
        for _col in columns:
            # Skips columns
            if _col in skip:
                continue
            # Include columns
            if _col in ['x', 'y', 'z']:
                _columns.append('ROUND(AVG({}), 8) AS `{}`'.format(_col, _col))
                if _col == 'x': _var_list.append('longitude')
                elif _col == 'y': _var_list.append('latitude')
                elif _col == 'z': _var_list.append('height')
            elif _col == 't':
                #FROM_UNIXTIME()
                _columns.append('AVG(UNIX_TIMESTAMP(t)) AS `t`')
                _var_list.append('time')
            elif _col == 'ms':
                _columns.append('MIN(ms) AS `ms`')
                _var_list.append('milliseconds')
            elif _col == 'quality':
                _columns.append('MAX(quality) as `quality`')
                _var_list.append('quality')
            else:
                _columns.append(_col)
                _var_list.append(_col)
            # Determines data type
            if _col in ['project_id', 'station_id']:
                _dtype.append('int')
            else:
                _dtype.append('float64')
        #
        _counter = 0
        if isinstance(var, basestring):
            var = [var]
        for _i, _item in enumerate(var):
            if isinstance(_item, list):
                _item, _var_suffix, _var_options = (_item[0], _item[1],
                    ' AND {}'.format(' AND '.join(_item[2:])))
            else:
                _var_suffix = _var_options = ''
            if isinstance(_item, basestring):
                _var = self._conn.datavars_read(name='{}'.format(_item),
                    debug=debug)
            elif isinstance(_item, atlantis.data.Variable):
                _var = _item
            else:
                raise ValueError(
                    'Invalid variable type ({})'.format(type(_item))
                )
            # Skips unknown or variables.
            if _var == None:
                continue
            if _counter > 0:
                _where += ' OR '
                _pivots += ', '
            _where += 'var_id={}'.format(_var.id)
            _pivots += 'AVG(CASE WHEN var_id={}{} THEN {} END) `{}{}`'.format(
                _var.id, _var_options, column, _var.standard_name, _var_suffix
            )
            _var_name = '{}{}'.format(_var.standard_name, _var_suffix)
            _var_list.append(_var_name)
            _var_return[_var_name] = _var
            _dtype.append('float64')
            _counter += 1
            # Add units to convert from and to.
            if units == 'units':
                if _var.units == None:
                    _var.units = _var.canonical_units
                _units[_var_name] = (_var.canonical_units, _var.units)
            elif units == 'canonical_units':
                _units[_var_name] = (_var.canonical_units,
                    _var.canonical_units)
            elif isinstance(units, list):
                _units[_var_name] = (_var.canonical_units, units[_i])
            elif isinstance(units, dict):
                _units[_var_name] = (_var.canonical_units,
                        units[_var.standard_name])

        # Appends level, project and station
        _where_project = ''
        if isinstance(project, Project):
            _where_project += 'project_id={}'.format(project.id)
        if isinstance(level, basestring):
            level = self._conn.datalevels_read(level, debug=debug)
        if isinstance(level, Level):
            if _where_project != '':
                _where_project = '{} AND '.format(_where_project)
            _where_project += 'level_id={}'.format(level.id)
        if isinstance(station, Station):
            if _where_project != '':
                _where_project = '{} AND '.format(_where_project)
            _where_project += '(project_id={} AND station_id={}) '.format(
                station.project.id, station.id)
        # Appends location `where` SQL clauses.
        _where = self._append_where(_where, 't', t, radius_t)
        _where = self._append_where(_where, 'z', z, radius_v)
        _where = self._append_where(_where, 'y', y, radius_h)
        _where = self._append_where(_where, 'x', x, radius_h)
        _where = self._append_where(_where, 't', None, tlim)
        _where = self._append_where(_where, 'z', None, zlim)
        _where = self._append_where(_where, 'y', None, ylim)
        _where = self._append_where(_where, 'x', None, xlim)
        if quality == 1:
            _where = _where + \
                ' AND (quality IS NULL OR quality=0 OR quality=1)'
        elif quality == 4:
            _where = _where + \
                ' AND (quality IS NULL OR quality=0 OR quality=1 OR quality=4)'
        elif quality == -1:
            _where = _where + \
                ' AND (quality=1)'
        elif quality == -4:
            _where = _where + \
                ' AND (quality=4)'
        else:
            raise ValueError('Quality flag `{}` not implemented yet.'.format(quality))
        # Merges where clauses
        if (where != None):
            if _where == '':
                _where = where
            else:
                _where = '({}) AND ({})'.format(where, _where)
        if (_where_project != ''):
            if _where == '':
                _where = _where_project
            else:
                _where = '({}) AND ({})'.format(_where_project, _where)

        # The almost final SQL command.
        _columns = ', '.join(_columns)
        _sql = 'SELECT {}, {} FROM datavalues WHERE ({}) GROUP BY {}'.format(
            _columns, _pivots, _where, groupby)
        # Includes data sorting
        if sort == True:
            _sql = '{} ORDER BY t, ms, x, y, z'.format(_sql)
        elif sort != '':
            _sql = '{} ORDER BY {}'.format(_sql, sort)
        # Fetches SQL query
        _query = self._conn.fetch(_sql, debug=debug)
        if len(_query) == 0:
            return None, None

        # Recast nested tuple to a python list and flatten it so it's a proper
        # iterable. Since we imported `unicode_literals` from `__future__`,
        # we have to change the data type of each item in `_var_list` to `str`.
        _var_list = [str(_item) for _item in _var_list]
        if masked:
            _query = ma.array(_query, dtype=zip(_var_list, _dtype))
            for _item in _var_list:
                _query[_item] = ma.masked_invalid(_query[_item].data)
        else:
            _query = array(_query, dtype=zip(_var_list, _dtype))
        # Converts UNIX timestamp to matplotlib date number.
        _query['time'] = epoch2num(_query['time'])

        # Converts the units of the data.
        for key, (_from, _to) in _units.items():
            _query[key] = atlantis.units.fromAtoB(_query[key].astype(float),
                _from, _to, vtype=column)
            _var_return[key].units = _to

        # And finally, returns the query results.
        return _query, _var_return


    def _asarray(self, x, size=None, dtype=None):
        """Converts any number of array like variable to numpy array."""
        try:
            if x is None:
                x = asarray([None])
            elif type(x) in [str, int, float, float32, float64]:
                x = asarray([x], dtype=type(x))
            elif type(x) in [list, tuple]:
                x = asarray(x)
            # Converts to numpy time, if appropriate
            if (dtype == 'time') & (type(x[0]) in [int, float, float32,
                float64]):
                x = asarray(num2date(x))
        except:
            print type(x), x, size, dtype
            raise ValueError('AAAhhhh!!!')
        if (size != None) & (x.size == 1):
            return x.repeat(size)
        else:
            return x


    def _read_asarray(self, x, i):
        """If x is an array, returns i-th value, otherwise returns x"""
        try:
            if isnan(x[i]) | is_masked(x[i]):
                return None
            else:
                return x[i]
        except:
            return x


    def _append_where(self, where, column, x, radius=0.):
        """Appends spatio-temporal where clause to where."""
        if where != '':
            if (where[0] != '(') & (where[-1] != ')'):
                where = '({})'.format(where)
        if isinstance(x, tuple) | isinstance(x, ndarray):
            raise ValueError('Invalid data type `{}`'.format(type(ndarray)))
        elif (column == 't'):
            _start, _stop = self._start_stop(x, radius)
            if type(_start) not in [list, tuple]:
                _start = [_start]
            if type(_stop) not in [list, tuple]:
                _stop = [_stop]
            where_time = ''
            for _istart, _istop in zip(_start, _stop):
                if where_time != '':
                    where_time += ' OR '
                if (_istart != None) & (_istop == None):
                    where_time += '({}={})'.format(column, _istart)
                elif _istop != None:
                    where_time += '({} BETWEEN \'{}\' AND \'{}\')'.format(
                        column, _istart, _istop
                    )
            if (where == '') & (where_time != ''):
                where = '({})'.format(where_time)
            elif where_time != '':
                where = '{} AND ({})'.format(where, where_time)
        elif type(x) in [int, float, str]:
            if (radius == 0) or (radius == None):
                where = '{} AND ({}={})'.format(where, column, x)
            else:
                where = '{} AND ({}>={}-{} AND {}<={}+{})'.format(
                    where, column, x, radius, column, x, radius
                )
        elif type(radius) in [list, tuple, ndarray]:
            where = '{} AND ({}>={} AND {}<={})'.format(
                where, column, radius[0], column, radius[-1]
            )
        #elif (x == None) and (radius == None):
        #    where = '{} AND ({} IS NULL)'.format(where, column)
        #
        return where


    def _start_stop(self, t, radius=0.):
        _fmt = '%Y-%m-%d %H:%M:%S'
        if t == None:
            if isinstance(radius, list) | isinstance(radius, tuple):
                if type(radius[0]) in [float64, float32, float, int, str,
                                       unicode]:
                    return (self._get_date(radius[0]).strftime(_fmt),
                        self._get_date(radius[-1]).strftime(_fmt))
                else:
                    return ([self._get_date(item[0]) for item in radius],
                        [self._get_date(item[1]) for item in radius])
            else:
                return None, None
        else:
            _date = self._get_date(t)
            if radius > 0.:
                _delta = timedelta(radius)
                return ((_date - _delta).strftime(_fmt),
                    (_date + _delta).strftime(_fmt))
            else:
                return _date.strftime(_fmt), None


    def _get_date(self, t):
        if type(t) in [int, float, float32, float64]:
            return num2date(t)
        elif isinstance(t, basestring):
            return parse(t)
        else:
            raise ValueError('Invalid data type.')
