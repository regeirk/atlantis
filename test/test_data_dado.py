# -*- coding: utf-8 -*-
"""Atlantis data framework.

Atlantis is a Python library for atmospheric, oceanographic and
hydrographic data analysis and visualization.

All analysis is centered around a common framework for structured data.
The package has to be able to handling multi-dimensional data and
associated metadata. Much of this is based uppon Iris library

This module runs tests on the data management module.

AUTHOR
    Sebastian Krieger
    email: sebastian.krieger@usp.br

REVISION
    1 (2014-09-12 17:55 -0300)

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

import unittest

from numpy import arange, array
from numpy.random import random

from atlantis.data import dado
from atlantis import data


class TestData(unittest.TestCase):
    def setUp(self):
        self.connection = dado.Connection(user='dado', password='CEBIMar',
            db='dado')
        self.sequence = dado.Sequence(self.connection)
        self.variable = data.get_standard_variable('sea_water_temperature')
        self.level_raw = dado.Level(id=0, name='raw', temporal='NULL',
            horizontal_x='NULL', horizontal_y='NULL', vertical='NULL',
            description='Raw (level 0) data.')
        self.project_dummy = dado.Project(id=0, parent=None, name='dummy',
            description='A dummy project for testing purposes.')
        self.station_dummy = dado.Station(project=self.project_dummy,
            id=0, name='free',
            description='A free station for the dummy project.')


    def tearDown(self):
        self.connection.close()
        self.connection = None


    def test_connection_execute(self):
        sql = 'SHOW TABLES;'
        self.connection.execute(sql)


    def test_connection_fetch(self):
        sql = 'SHOW TABLES;'
        a = self.connection.fetch(sql)
        self.assertEqual('datalevels', a[0][0])


    def test_connection_datavars_insert(self):
        var_list = ['sea_water_temperature', 'sea_water_salinity', 'depth']
        for var in var_list:
            self.connection.datavars_insert(data.get_standard_variable(var))


    def test_connection_datavars_read(self):
        a = self.connection.datavars_read(where='var_id=21001')
        b = self.connection.datavars_read(where=['var_id=21001',
            'name="sea_water_temperature"'])
        self.assertEqual('sea_water_temperature', a.standard_name)
        self.assertEqual('sea_water_temperature', b.standard_name)


    def test_connection_datalevels_insert(self):
        self.connection.datalevels_insert(self.level_raw)


    def test_connection_datalevels_read(self):
        a = self.connection.datalevels_read(where='level_id=0')
        self.assertEqual('raw', a.name)


    def test_connection_project_insert(self):
        self.connection.projects_insert(self.project_dummy)


    def test_connection_project_read(self):
        a = self.connection.projects_read(where='project_id=0')
        self.assertEqual('dummy', a.name)


    def test_connection_station_insert(self):
        self.connection.stations_insert(self.station_dummy)


    def test_connection_station_read(self):
        a = self.connection.stations_read(where='project_id=0')
        self.assertEqual('free', a.name)
        self.assertEqual(0, a.project.id)


    def test_sequence_write_onevar_oneentry(self):
        station = self.station_dummy
        level = self.level_raw
        var = self.connection.datavars_read(
            where='name="{:s}"'.format(self.variable.standard_name)
        )
        # Location of BATS site.
        x, y, z = 31.666666666666668, -64.16666666666667, -3.
        t = '1978-06-21 17:44'
        ms = 500
        val = random()
        #
        self.sequence.write(t, z, y, x, val, station, level, var, ms=ms)


    def test_sequence_write_onevar_multipleentries(self):
        station = self.station_dummy
        level = self.level_raw
        var = self.connection.datavars_read(
            where='name="{:s}"'.format(self.variable.standard_name)
        )
        # Location of BATS site.
        x, y, z = 31.666666666666668, -64.16666666666667, -3.
        t0 = 722256.7388888889 # '1978-06-21 17:44'
        t = arange(25) + t0
        ms = 500
        val = random(25)
        #
        self.sequence.write(t, z, y, x, val, station, level, var, ms=ms)


    def test_sequence_write_multiplevars_multipleentries(self):
        station = self.station_dummy
        level = self.level_raw
        var = self.connection.datavars_read()
        # Location of BATS site.
        x, y, z = 31.666666666666668, -64.16666666666667, -3.
        t0 = 722256.7388888889 # '1978-06-21 17:44'
        t = arange(10) + t0
        ms = 500
        val = []
        for v in var:
            val.append(array(random(10)))
        #
        self.sequence.write(t, z, y, x, val, station, level, var, ms=ms)


    def test_sequence_read(self):
        dat = self.sequence.read(['sea_water_temperature',
            'sea_water_salinity', 'pressure'])


    def test_sequence_read_near_BATS(self):
        t = 722256.7388888889
        x, y = 31.666666666666668, -64.16666666666667
        dat = self.sequence.read('sea_water_temperature', t=t, y=y, x=x,
            radius_h=1., radius_t=9)


    #def test_dataset_add_variable(self):
    #    self.dataset.add_variable(self.variable)


    #def test_dataset_get_variable(self):
    #    self.dataset.add_variable(self.variable)
    #    a = self.dataset.get_variable()
    #    self.assertEqual(self.variable.standard_name,
    #        a[self.variable.standard_name].standard_name)

    #def test_connection_dataset_insert(self):
    #    # First we have to create a proper dataset
    #    project = self.connection.projects_read(where='project_id=0')
    #    level = self.connection.datalevels_read(where='level_id=0')
    #    ds = dado.Dataset(project=project, level=level, name='test',
    #        description='This is a test dataset for the dummy project.')
    #    var_list = ['sea_water_temperature', 'sea_water_salinity', 'depth']
    #    for var in var_list:
    #        self.connection.datavars_insert(data.get_standard_variable(var))
    #    # And now we add it to the database.
    #    self.connection.dataset_insert(self, dataset)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
