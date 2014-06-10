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
    2 (2013-10-26 18:47 -0300 DST)
    1 (2013-06-27 18:56 -0300)

"""
from __future__ import division

__version__ = '$Revision: 1 $'
# $Source$

import unittest
import numpy

try:
    reload(data)
except:
    from atlantis import data

import klib

class TestData(unittest.TestCase):
    def setUp(self):
        dx = dy = 0.25
        T0 = 693596.5 # '1900-01-01 12:00'
        self.grid_time = numpy.arange(10) + T0
        self.grid_height = numpy.arange(0)
        self.grid_lat = numpy.arange(-90., 90., dy) + dy / 2.
        self.grid_lon = numpy.arange(0., 360., dx) + dx / 2.
        self.grid_N = numpy.arange(self.grid_time.size)
        self.grid_K = numpy.arange(self.grid_height.size)
        self.grid_J = numpy.arange(self.grid_lat.size)
        self.grid_I = numpy.arange(self.grid_lon.size)
        #
        x, y = self.grid_lon, self.grid_lat
        xx, yy = numpy.meshgrid(x, y)
        t = 0
        k = 2 * numpy.pi / 180
        l = 2 * numpy.pi / 720
        w = 2 * numpy.pi / 365.25
        self.grid_SST = numpy.cos(k * xx * l * yy - w * t)
        self.grid_dSSTdx = - k * l * yy * numpy.sin(k * xx * l * yy - w * t)
        self.grid_d2SSTdx2 = -(k*l*yy)**2 * numpy.cos(k * xx * l * yy - w * t)
        self.grid_dSSTdy = - k * l * xx * numpy.sin(k * xx * l * yy - w * t)
        self.grid_d2SSTdy2 = -(k*l*xx)**2 * numpy.cos(k * xx * l * yy - w * t)
        self.grid_grad = self.grid_dSSTdx + 1j * self.grid_dSSTdy
        #
        l = 2 * numpy.pi / 45
        self.grid_UV = (numpy.sin(k * xx) + numpy.cos(l * yy) +
            1j * (numpy.cos(k * xx) + numpy.sin(l * yy)))
        self.grid_curl = -k * numpy.sin(k * xx) - (- l * numpy.sin(l * yy))
        self.grid_div = k * numpy.cos(k * xx) + l * numpy.cos(l * yy)
        #
        self.grid = data.Grid()
        self.grid.name = 'Test dataset'
        self.grid.dimensions = dict(
            n = self.grid_time.size,
            k = self.grid_height.size,
            j = self.grid_lat.size,
            i = self.grid_lon.size
        )
        var_list = ['time', 'height', 'latitude', 'longitude',
            'sea_surface_temperature']
        for var in var_list:
            self.grid.variables[var] = data.get_standard_variable(var)
        self.grid.variables['time'].data = self.grid_time
        self.grid.variables['time'].canonical_units = ('days since 0001-01-01'
            ' UTC')
        self.grid.variables['height'].data = self.grid_height
        self.grid.variables['latitude'].data = self.grid_lat
        self.grid.variables['longitude'].data = self.grid_lon
        self.grid.variables['sea_surface_temperature'] = self.grid_SST


    def tearDown(self):
        self.grid = None
        self.grid_time = None
        self.grid_height = None
        self.grid_lat = None
        self.grid_lon= None
        self.grid_N = None
        self.grid_K = None
        self.grid_J = None
        self.grid_I = None

    
    def test_get_standard_variable(self):
        var = data.get_standard_variable('sea_surface_temperature')
        self.assertEqual('sea_surface_temperature', var.standard_name)
        self.assertEqual('K', var.canonical_units)


    def test_get_standard_variable_invalid(self):
        # Should raise an exception for invalid standard_name
        self.assertRaises(ValueError, data.get_standard_variable,
            'sea_surface_height')


    def test_grid_curl(self):
        #
        # Calculates the curl of a vector map
        curl = self.grid.curl(self.grid_UV, q=9, units='deg', cyclic=False)
        #
        # Plots the results
        fig = klib.graphics.figure(fp=dict(), orientation='worldmap')
        klib.gis.map(self.grid_lon, self.grid_lat, curl, z2=self.grid_UV,
            da=[51, 51], fig=fig, show=True, profile=False)
        fig.savefig('./test_data_grid_curl.png', dpi=150)
        #
        # Do the tests
        numpy.testing.assert_array_almost_equal(curl, self.grid_curl,
            decimal=4)


    def test_grid_div(self):
        #
        # Calculates the curl of a vector map
        div = self.grid.div(self.grid_UV, q=9, units='deg', cyclic=False)
        #
        # Plots the results
        fig = klib.graphics.figure(fp=dict(), orientation='worldmap')
        klib.gis.map(self.grid_lon, self.grid_lat, div, z2=self.grid_UV,
            da=[51, 51], fig=fig, show=True, profile=False)
        fig.savefig('./test_data_grid_divergence.png', dpi=150)

        #
        # Do the tests
        numpy.testing.assert_array_almost_equal(div, self.grid_div,
            decimal=4)


    def test_grid_grad(self):
        #
        # Calculates the gradient of a scalar map
        grad = self.grid.grad(self.grid_SST, q=9, units='deg', cyclic=False)
        #
        # Plots the results
        fig = klib.graphics.figure(fp=dict(), orientation='worldmap')
        klib.gis.map(self.grid_lon, self.grid_lat, self.grid_SST, z2=grad,
            da=[51, 51], fig=fig, show=True, profile=False)
        fig.savefig('./test_data_grid_gradient.png', dpi=150)
        
        #
        # Do the tests
        numpy.testing.assert_array_almost_equal(grad, self.grid_grad,
            decimal=4)
    
    
    def test_grid_derivative_x(self):
        #
        # Calculates first and second derivative
        dSSTdx = self.grid.derivative(self.grid_SST, axis='longitude', q=9,
            p=1, cyclic=False)
        d2SSTdx2 = self.grid.derivative(self.grid_SST, axis='longitude', q=9,
            p=2, cyclic=False)
        #
        # Plots the results
        ap = dict(left=0.1, bottom=0.12, right=0.82, top=0.95, wspace=0.1,
            hspace=0.5)
        fig = klib.graphics.figure(fp=dict(), ap=ap, orientation='landscape')
        crange0 = numpy.arange(-1, 1.1, .1)
        klib.graphics.contour(self.grid_lon, self.grid_lat, self.grid_SST,
            fig=fig, subplot=(3, 1, 1), orientation='vertical',
            crange=crange0, label='a) $f(x, y)$')
        RMSD = klib.common.latex_scientific(((dSSTdx.data -
            self.grid_dSSTdx)**2).mean())
        crange1 = 0.1 * numpy.arange(-1, 1.1, .1)
        ax = klib.graphics.contour(self.grid_lon, self.grid_lat, dSSTdx.data,
            fig=fig, subplot=(3, 1, 2), orientation='vertical', crange=crange1,
            scale=1e-1,
            label=r'b) $\frac{\partial f}{\partial x}$, $RMSD=%s$' % (RMSD))
        ax.contour(self.grid_lon, self.grid_lat, self.grid_dSSTdx,
            crange1[::5], colors='k')
        RMSD = klib.common.latex_scientific(((d2SSTdx2.data -
            self.grid_d2SSTdx2)**2).mean())
        crange2 = 0.01 * numpy.arange(-1, 1.1, .1)
        bx = klib.graphics.contour(self.grid_lon, self.grid_lat, d2SSTdx2.data,
            fig=fig, subplot=(3, 1, 3), orientation='vertical', crange=crange2,
            scale=1e-2,
            label=r'c) $\frac{\partial^2 f}{\partial x^2}$, $RMSD=%s$' %
            (RMSD))
        bx.contour(self.grid_lon, self.grid_lat, self.grid_d2SSTdx2,
            crange2[::5], colors='k')
        fig.savefig('./test_data_grid_derivative_x.png', dpi=150)
        #
        # Do the tests
        numpy.testing.assert_array_almost_equal(dSSTdx, self.grid_dSSTdx,
            decimal=4)
        numpy.testing.assert_array_almost_equal(d2SSTdx2, self.grid_d2SSTdx2,
            decimal=4)


    def test_grid_derivative_y(self):
        #
        # Calculates first and second derivative
        dSSTdy = self.grid.derivative(self.grid_SST, axis='latitude', q=9,
            p=1, cyclic=False)
        d2SSTdy2 = self.grid.derivative(self.grid_SST, axis='latitude', q=9,
            p=2, cyclic=False)
        #
        # Plots the results
        ap = dict(left=0.1, bottom=0.12, right=0.82, top=0.95, wspace=0.1,
            hspace=0.5)
        fig = klib.graphics.figure(fp=dict(), ap=ap, orientation='landscape')
        crange0 = numpy.arange(-1, 1.1, .1)
        klib.graphics.contour(self.grid_lon, self.grid_lat, self.grid_SST,
            fig=fig, subplot=(3, 1, 1), orientation='vertical',
            crange=crange0, label='a) $f(x, y)$')
        RMSD = klib.common.latex_scientific(((dSSTdy.data -
            self.grid_dSSTdy)**2).mean())
        crange1 = 0.1 * numpy.arange(-1, 1.1, .1)
        ax = klib.graphics.contour(self.grid_lon, self.grid_lat, dSSTdy.data,
            fig=fig, subplot=(3, 1, 2), orientation='vertical', crange=crange1,
            scale=1e-1,
            label=r'b) $\frac{\partial f}{\partial y}$, $RMSD=%s$' % (RMSD))
        ax.contour(self.grid_lon, self.grid_lat, self.grid_dSSTdy,
            crange1[::5], colors='k')
        RMSD = klib.common.latex_scientific(((d2SSTdy2.data -
            self.grid_d2SSTdy2)**2).mean())
        crange2 = 0.01 * numpy.arange(-1, 1.1, .1)
        bx = klib.graphics.contour(self.grid_lon, self.grid_lat, d2SSTdy2.data,
            fig=fig, subplot=(3, 1, 3), orientation='vertical', crange=crange2,
            scale=1e-2,
            label=r'c) $\frac{\partial^2 f}{\partial y^2}$, $RMSD=%s$' %
            (RMSD))
        bx.contour(self.grid_lon, self.grid_lat, self.grid_d2SSTdy2,
            crange2[::5], colors='k')
        fig.savefig('./test_data_grid_derivative_y.png', dpi=150)
        #
        # Do the tests
        numpy.testing.assert_array_almost_equal(dSSTdy, self.grid_dSSTdy,
            decimal=4)
        numpy.testing.assert_array_almost_equal(d2SSTdy2, self.grid_d2SSTdy2,
            decimal=4)

def main():
    klib.graphics.pylab.close('all')
    unittest.main()


if __name__ == '__main__':
    main()
