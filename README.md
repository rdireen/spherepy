SpherePy
========

SpherePy is a package for working with scalar and vector spherical harmonics.
It provides:

* scalar and vector spherical harmonic decompositions
* objects for algebraically manipulating coefficients
* the ability to plot coefficients and patterns

Badges
------

**AppVeyor:**

[![Build status](https://ci.appveyor.com/api/projects/status/ccwuv424wao4rbly?svg=true)](https://ci.appveyor.com/project/rdireen/spherepy)

**TravisCI:**

[![Build Status](https://travis-ci.org/rdireen/spherepy.svg?branch=master)](https://travis-ci.org/rdireen/spherepy)
[![Coverage Status](https://coveralls.io/repos/rdireen/spherepy/badge.svg?branch=master)](https://coveralls.io/r/rdireen/spherepy?branch=master)

**PyPI:**

[![Latest Version](https://pypip.in/version/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)
[![Supported Python versions](https://pypip.in/py_versions/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)
[![Supported Python implementations](https://pypip.in/implementation/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)
[![Development Status](https://pypip.in/status/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)
[![Download format](https://pypip.in/format/spherepy/badge.svg)](https://pypi.python.org/pypi/spherepy/)



Better Documentation
--------------------

I'm working on a more complete set of documentation 
[HERE](http://www.direentech.com/docs/spherepy).

	
Install
=======


**Ubuntu:**
Before installing SpherePy you must install build-essential, python-dev, and NumPy

```
$ sudo apt-get install build-essential python-dev
```

For NumPy you need to decide if you want to build it yourself with

```
$ sudo pip install numpy
```

or download the package with

```
$ sudo apt-get install python-numpy
```

I have been building NumPy using the pip method, but it takes a long time. 

Then you can 

```
$ sudo pip install spherepy
```
    
**Windows:**
Make sure you have Numpy on your machine, then

```
$ pip install spherepy
```
	
Plotting
========

If you would like to use the plotting routines within SpherePy, install matplotlib:

```
$ sudo pip install matplotlib
```
	
Quick Example
=============

```python
import spherepy as sp 
c = sp.random_coefs(4, 4) # generate some random coefficients
sp.pretty_coefs(c)
p = sp.ispht(c, 50, 50) # inverse spherical transform to pattern
sp.plot_sphere_mag(p) # plot the magnitude of the pattern
```

Credits
=======
The algorithms within this package have been implemented, in part, using the documentation within 
the NIST Interagency/Internal Report (NISTIR) - 3955 [citation](http://www.nist.gov/manuscript-publication-search.cfm?pub_id=1051).
The majority or the code, however, has been developed using Ronald C. Wittmann's unpublished notes.

Contributing
============
Reporting bugs, suggesting features, helping with documentation, and adding to the code is very welcome. See
[Contributing](CONTRIBUTING.md). 

License
=======

Copyright (C) 2015  Randy Direen <spherepy@direentech.com>.
SpherePy is licensed under GNU General Public License, version 3, a copy of this license has been provided within the COPYING file in this directory, and can also be found at <http://www.gnu.org/licenses/>.
 
