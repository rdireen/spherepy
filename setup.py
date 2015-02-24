#!/usr/bin/env python

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
          'Programming Language :: Python :: Implementation :: CPython'
          'Topic :: Scientific/Engineering :: Mathematics'
      ],
      url='https://github.com/rdireen/spherepy',  # url to github repo
      download_url='https://github.com/rdireen/spherepy/tarball/0.1',
      license='GPLv3',
      install_requires=['numpy', 'setuptools', 'six'],
      packages=['spherepy'],
      package_dir={'spherepy':'spherepy', 'test':'spherepy/test'},
      package_data={'spherepy':['pkg_info.json']},
      include_package_data=True,
      test_suite='nose.collector',
      tests_require=['nose'],
      ext_modules=[csphi_module]  
     )

