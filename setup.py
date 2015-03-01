#!/usr/bin/env python

# Copyright (C) 2015  Randy Direen <spherepy@direentech.com>
#
# This file is part of SpherePy.
#
# SpherePy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SpherePy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SpherePy.  If not, see <http://www.gnu.org/licenses/>

"""

Note about NumPy:

    NumPy is a necessary part of SpherePy, but there are many ways to install
    it. On Ubuntu you can install the the python-numpy package with apt-get or
    use pip, and on Windows you can use full systems like Anaconda, download
    wheels from unofficial sites, or if you have the compilers for the right 
    version of the Python distribution you can build it yourself. Since there
    are so many ways to get NumPy, I think I will require people to install
    and update it themselves and not put it in the requirements here.

"""

import os
import sys
import json

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'spherepy/pkg_info.json')) as fp:
    _info = json.load(fp)

def readme():
    with open('README.md') as f:
        return f.read()

__version__ = _info['version']
__author__ = _info['author']
__email__ = _info['email']

try:
    from setuptools import setup, Extension
except ImportError:
    print("SpherePy requires setuptools in order to build. Install " + \
          "setuptools using your package manager (possibly " + \
          "python-setuptools) or using pip (i.e., pip install "
          "setuptools")
    sys.exit(1)
 
try:   
    import numpy
    # Obtain the numpy include directory.  
    # This logic works across numpy versions.
    try:
        numpy_include = numpy.get_include()
    except AttributeError:
        numpy_include = numpy.get_numpy_include()
    
except ImportError:
    print("SpherePy requires NumPy for compiling c extensions. Install " + \
          "NumPy using your packag manager (possibly, python-numpy) or " + \
          "using pip (i.e., pip install numpy).")
    sys.exit(1)

description = 'Numerical routines for working with spherical harmonic ' + \
              'coefficients' 
 
srcs = ['src/csphi.c', 'src/csphi_wrap.c', 'src/kiss_fft.c']  
headers=['src/csphi.h', 'src/kiss_fft.h', 'src/_kiss_fft_guts.h']          
csphi_module = Extension('_csphi',
                         sources=srcs,
                         include_dirs=['src', numpy_include]
                       )

setup(name='spherepy',
      version=__version__,
      author=__author__,
      author_email=__email__,
      description=description,
      long_description=readme(),
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Operating System :: Microsoft :: Windows :: Windows 7',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: Implementation :: CPython',
          'Topic :: Scientific/Engineering :: Mathematics'
      ],
      url='https://github.com/rdireen/spherepy',  # url to github repo
      download_url='https://github.com/rdireen/spherepy/tarball/0.2',
      license='GPLv3',
      install_requires=['setuptools', 'six'],
      keywords=['sphere transform'],
      packages=['spherepy'],
      package_dir={'spherepy':'spherepy', 'test':'spherepy/test'},
      package_data={'spherepy':['pkg_info.json']},
      include_package_data=True,
      test_suite='nose.collector',
      tests_require=['nose'],
      ext_modules=[csphi_module],
      headers=headers  
     )

