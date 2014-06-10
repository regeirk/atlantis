import unittest
import numpy

class TestDataDummy(unittest.TestCase):
    def setUp(self):
        self.A = numpy.array([0, 1, 2, 3, 4, 5])
        self.B = numpy.array([0, 1, 2, 3, 4, 5])

    def test_AB(self):
        numpy.testing.assert_array_almost_equal(self.A, self.B)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
