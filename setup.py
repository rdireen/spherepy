#!/usr/bin/env python

import os
import sys
from glob import glob
import json

here = os.path.abspath(os.path.dirname(__file__))

with open('spherepy/pkg_info.json') as fp:
    _info = json.load(fp)

def readme():
    with open('README.md') as f:
        return f.read()

__version__ = _info['version']
__author__ = _info['author']

try:
    from setuptools import setup, find_packages
except ImportError:
    print("SpherePy requires setuptools in order to build. Install " + \
          "setuptools using your package manager (possibly " + \
          "python-setuptools) or using pip (e.i., pip install "
          "setuptools")
    sys.exit(1)

description = 'Numerical routines for working with spherical harmonic ' + \
              'coefficients' 

setup(name='spherepy',
      version=__version__,
      author="Randy Direen",
      author_email='spherepy.direentech.com',
      description=description,
      long_description=readme(),
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Mathematics'
      ],
      url='http://github.com/rdireen/spherepy',
      license='GPLv3',
      install_requires=['numpy','setuptools','matplotlib'],
      packages=['spherepy'],
      package_dir={'spherepy':'spherepy','test':'spherepy/test'},
      package_data={'spherepy':['pkg_info.json']},
      include_package_data = True,
      test_suite='nose.collector',
      tests_require=['nose']   
     )

