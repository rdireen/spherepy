.. spherepy documentation master file, created by
   sphinx-quickstart on Sat Feb  7 21:35:42 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
 
 
.. COMMENT: Put all references to media here
 

  

**********************


.. image:: images/logo_spherepy.png
   :height: 200px
   :align: right

   
SpherePy
********

.. warning:: 

   In mega alpha mode. Please fork this on GitHub (`spherepy <http://github.com/rdireen/spherepy>`_)
   and help me. 


**Jump straight to** :doc:`GettingStarted` **or the** :ref:`ref-toc`


SpherePy does spherical harmonic transforms. Given a
function on a sphere :math:`f(\theta, \phi)`, SpherePy breaks the function into a finite 
sum of spherical harmonics :math:`Y_{nm}(\theta, \phi)`:

.. math::
    f(\theta, \phi) = \sum_{n=0}^{N} \sum_{m=-n}^{n} c_{nm} Y_{nm}(\theta, \phi)

The output of the transform is the set of coefficients :math:`c_{nm}`. 

Here is a simple code example::

    >>> import spherepy as sp
    >>> f = sp.random_patt_uniform(50,50) # create a random pattern
    >>> c = sp.spht(f, 20, 20)  # spherical harmonic transform of f
    >>> sp.pcolor_coefs(c) # plot the set of coefficients (see below)
    >>> for n in range(0, 2): # print out first couple of modes
    ...     for m in range(-n, n + 1):
    ...         print("c[{0}, {1}] = {2}".format(n, m, c[n, m])) 
    ...
    c[0, 0] = (1.2-0.5j)
    c[1, -1] = (-1.3+3.0j)
    c[1, 0] = (0.1-0.8j)
    c[1, 1] = (2.3-2.2j) 
    >>> f2 = sp.ispht(c[0:4,:], 50, 50) # do the inverse transform of the first 5 modes
    >>> sp.plot_sphere_mag(f2) # plot the pattern (see below)

and here are the plots:

.. image:: images/rnd_coefsp.png
   :width: 400px
   
.. image:: images/rnd1p.png
   :width: 350px


Now What?
=========

**Get SperePy**: :doc:`GettingStarted`

.. _ref-toc:

Table of Contents
=================
.. figure:: images/me.jpg
   :width: 250px
   :align: right
   :figclass: align-center
   
   A picture of me

.. toctree::
   :maxdepth: 2
   
   GettingStarted
   Tutorial1
   GruesomeDetails
   Tutorial2
   Details
   CodeDoc
   
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`









   


   





