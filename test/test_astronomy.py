import unittest
import numpy

from atlantis import astronomy

class TestAstronomy(unittest.TestCase):
    def setUp(self):
        self.K = astronomy.constants()
        #
        self.lon_1 = numpy.zeros(10)
        self.lat_1 = numpy.arange(-5, 5)
        self.x_1, self.y_1 = numpy.meshgrid(self.lon_1, self.lat_1 * self.K.b)

    def test_metergrid(self):
        x, y = astronomy.metergrid(self.lon_1, self.lat_1, unit='m')
        numpy.testing.assert_array_almost_equal(self.x_1, x)
        numpy.testing.assert_array_almost_equal(self.y_1, y)

    def test_metergrid_km(self):
        x, y = astronomy.metergrid(self.lon_1, self.lat_1, unit='km')
        numpy.testing.assert_array_almost_equal(self.x_1 * 1e-3, x)
        numpy.testing.assert_array_almost_equal(self.y_1 * 1e-3, y)

    def test_metergrid_nm(self):
        x, y = astronomy.metergrid(self.lon_1, self.lat_1, unit='nm')
        numpy.testing.assert_array_almost_equal(self.x_1 / 1852, x)
        numpy.testing.assert_array_almost_equal(self.y_1 / 1852, y)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
