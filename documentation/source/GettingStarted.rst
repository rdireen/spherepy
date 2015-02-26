.. include global.rst

Getting Started
***************

Works on Python 2.7, 3.3, and 3.4 under Linux and Windows (not sure about OS X).



Install
=======

* Setting up :ref:`ref-pip`
* Install :ref:`ref-linux`
* Install :ref:`ref-windows`
* Install :ref:`ref-source`
* :ref:`ref-uninstall`

.. _ref-pip:

pip
---
Whether you are on Linux, OS X, or Windows, make sure you have the latest
version of `pip <https://pip.pypa.io/en/latest/installing.html>`_ installed 
(https://pip.pypa.io/en/latest/installing.html).
 
**Upgrading**:

On Linux or OS X::

    $ pip install -U pip

On Windows::

    $ python -m pip install -U pip



.. _ref-linux:

SpherePy on Linux
-----------------

SpherePy requires C extensions to be compiled. On Ubuntu, you will have to do the following::

    $ sudo apt-get install build-essentials python-dev

Building C extensions also requires NumPy to be installed::

    $ sudo pip install numpy

You may then install SpherePy::

    $ sudo pip install spherepy

You will most likely also want matplotlib::

    $ sudo pip install matplotlib
	
	
	
	
.. _ref-windows:

SpherePy on Windows
-------------------

SpherePy requires NumPy, which contains a lot of binaries. There are a lot of ways to get 
NumPy installed on Windows, but one of the easiest ways is to install a Python distribution
such as `Anaconda <http://continuum.io/downloads>`_. 

Now we should be able to install SpherePy::	

    $ pip install spherepy
	
We can also install matplotlib to plot things::

    $ sudo pip install matplotlib


.. _ref-source:

From Source Code
----------------

The source code is hosted at GitHub. For the latest (but not necessarily most 
stable) version of the code use git::

    $ git clone https://github.com/rdireen/spherepy.git

You can visit the project directly at https://github.com/rdireen/spherepy.
Fork my code there and start helping out; I don't really know what I'm doing.

To install from source, you will need tools for building the python C extensions. 
There is a version of the code completely implemented in Python, but I couldn't get
the low level routines to run fast. To compile C extensions on 
Ubuntu you'll need build-essentials and python-dev::

    $ apt-get install build-essentials python-dev

.. note::
    I've been able to develop these extensions on Windows as well by using 
    the compiler `here <http://www.microsoft.com/en-us/download/details.aspx?id=44266>`_.

Building C extensions also requires NumPy to be installed::

    $ pip install numpy

You should now be able to install the code by entering the directory containing the 
*setup.py* script and typing::

    $ python setup.py install

To see if things are working, start python and type::

    >>> import spherepy



To plot stuff you'll need matplotlib::

    $ pip install matplotlib

	
.. _ref-uninstall:

Uninstall
---------

If you don't like it, you can cleanly remove SpherePy from your machine with::

    $ pip uninstall spherepy



Quick Example
=============

Plotting individual scalar spherical harmonics::

    >>> import spherepy as sp
    >>> C = sp.zeros_coefs(5,5)
    >>> C[2,0] = 1
    >>> p = sp.ispht(C, 50, 50)
    >>> sp.plot_sphere_mag(p)

As a result you should see a plot of the spherical function :math:`Y_{2,0}(\theta, \phi)`:

.. figure::  images/sph2_0p.png 
   :width: 600px





