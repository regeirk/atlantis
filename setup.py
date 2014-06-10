#!/usr/bin/env python

from distutils.core import setup

setup(
    name='atlantis',
    version='0.1',
    description=('Python library for atmospheric, oceanographic and '
        'hydrographic data analysis and visualization'),
    author='Sebastian Krieger',
    author_email='solutions@nublia.com',
    url='http://nublia.com/solutions/atlantis/',
    packages=['data', 'dynamics', 'graphics', 'signal', 'stats', 'tests'],
)
