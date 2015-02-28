SpherePy
========

SpherePy is a package for working with scalar and vector spherical harmonics.
It provides:

	* scalar and vector spherical harmonic decompositions
	* objects for algebraically manipulating harmonic coefficients
	* the ability to plot coefficients and patterns

**Pin Crazy:** 
	
[![Build status](https://ci.appveyor.com/api/projects/status/ccwuv424wao4rbly?svg=true)](https://ci.appveyor.com/project/rdireen/spherepy)
[![Build Status](https://travis-ci.org/rdireen/spherepy.svg?branch=master)](https://travis-ci.org/rdireen/spherepy)
[![Coverage Status](https://coveralls.io/repos/rdireen/spherepy/badge.svg?branch=master)](https://coveralls.io/r/rdireen/spherepy?branch=master)

[![Downloads](https://pypip.in/download/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)
[![Latest Version](https://pypip.in/version/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)
[![Supported Python versions](https://pypip.in/py_versions/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)
[![Supported Python implementations](https://pypip.in/implementation/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)

[![Development Status](https://pypip.in/status/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)
[![Wheel Status](https://pypip.in/wheel/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)
[![Egg Status](https://pypip.in/egg/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)

[![Download format](https://pypip.in/format/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)
[![License](https://pypip.in/license/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)



	
Install
=======

**pip**

For Windows I use wheels to distribute the binaries, so make sure pip is up to date.

**Ubuntu:**
Before installing SpherePy you must install build-essential, python-dev, and NumPy

    $ sudo apt-get install build-essential python-dev
    $ sudo pip install numpy

Then you can 

    $ sudo pip install spherepy
    
**Windows:**
Make sure you have Numpy on your machine, then

    $ pip install spherepy
	
Plotting
========

If you would like to use the plotting routines within SpherePy, install matplotlib:

    $ sudo pip install matplotlib
	
Quick Example
=============

    >>> import spherepy as sp 
    >>> c = sp.random_coefs(4, 4) # generate some random coefficients
    >>> p = sp.ispht(c, 50, 50) # inverse spherical transform to pattern
    >>> sp.plot_sphere_mag(p) # plot the magnitude of the pattern

License
=======

Copyright (C) 2015  Randy Direen <spherepy@direentech.com>.
SpherePy is licensed under GNU General Public License, version 3, a copy of this license has been provided within the COPYING file in this directory, and can also be found at <http://www.gnu.org/licenses/>.
 
