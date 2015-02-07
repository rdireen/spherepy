#!/usr/bin/env python

import os
import sys
import json

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here,'spherepy/pkg_info.json')) as fp:
    _info = json.load(fp)

def readme():
    with open('README.md') as f:
        return f.read()

__version__ = _info['version']
__author__ = _info['author']
__email__ = _info['email']

try:
    from setuptools import setup
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
      author=__author__,
      author_email=__email__,
      description=description,
      long_description=readme(),
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Mathematics'
      ],
      url='https://github.com/rdireen/spherepy', #url to github repo
      download_url = 'https://github.com/rdireen/spherepy/tarball/0.1',
      license='GPLv3',
      install_requires=['numpy', 'setuptools'],
      packages=['spherepy'],
      package_dir={'spherepy':'spherepy','test':'spherepy/test'},
      package_data={'spherepy':['pkg_info.json']},
      include_package_data = True,
      test_suite='nose.collector',
      tests_require=['nose']   
     )

