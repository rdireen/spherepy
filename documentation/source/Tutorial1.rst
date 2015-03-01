Working with Spherical Harmonics
********************************

.. warning:: 

   In mega alpha mode. Please fork this on GitHub (`spherepy <http://github.com/rdireen/spherepy>`_)
   and help me. 


Spherical harmonics :math:`Y_{nm}(\theta, \phi)` are defined for nonnegative :math:`n` and 
:math:`-n \leq m \leq n`. To work with these, you can create a set of coefficients :math:`c_{nm}`, where 
each element is set to zero::

    >>> import spherepy as sp
    >>> c = sp.zeros_coefs(2, 2)
    >>> sp.pretty_coefs(c)
    
    c[n,m]
    =====

    2:       0j             0j             0j             0j             0j 
    1:                      0j             0j             0j  
    0:                                     0j    
    n  -------------  -------------  -------------  -------------  -------------  
           m = -2         m = -1         m = 0          m = 1          m = 2    

.. sidebar:: Documentation

   Here is the documentation for the functions used above: :ref:`fun-zeroscoefs`,
   :ref:`fun-pcoefs`. If you execute :samp:`type(c)`, you'll notice that :samp:`c` 
   is a :ref:`fun-scalarcoefs`
   
The object :samp:`c` has 9 elements, the largest :samp:`n` value (referred to as :samp:`nmax`) is 2 and the largest 
:samp:`m` value (referred to as :samp:`mmax`) is 2.

|



|
Individual elements of :samp:`c` can be set if indexed properly::

    >>> c[1, -1] = 1.0 
    >>> c[2, 1] = -1.0j
    >>> sp.pretty_coefs(c)
    
    c[n,m]
    =====

    2:       0j             0j             0j            -1j             0j 
    1:                     1+0j            0j             0j  
    0:                                     0j    
    n  -------------  -------------  -------------  -------------  -------------  
           m = -2         m = -1         m = 0          m = 1          m = 2   
    
In the above example we have :math:`c_{1,-1}=1` and :math:`c_{2,1}=1j`. If we want to calculate the pattern 
:math:`f(\theta, \phi)`, 

.. math::
    f(\theta, \phi) = \sum_{n=0}^{2} \sum_{m=-n}^{n} c_{nm} Y_{nm}(\theta, \phi)

we call the inverse spherical harmonic transform routine :ref:`fun-ispht` ::

    >>> f = sp.ispht(c, 50, 80)
    >>> sp.plot_sphere_mag(f)

.. sidebar:: Documentation

   The transform is documented here :ref:`fun-ispht`, and note that :samp:`type(f)` is 
   a :ref:`fun-scalarpatt`. 

.. image:: images/sph_small_sum.png
   :width: 350px
   :align: center
   
The pattern :samp:`f` is a numerical representation of :math:`f(\theta, \phi)`. The NumPy array contained within :samp:`f` 
is a complex valued array with 50 rows and 80 columns::

    >>> f.nrows
    50
    >>> f.ncols
    80
    >>> narray = f.cdata # cdata is a NumPy array of type complex128
    
.. note::
   It might seem funny that :ref:`fun-ispht` returns a :ref:`fun-scalarpatt` object rather than a simple NumPy array. 
   I do things this way for two reason: first, the pattern itself needs to be preprocessed to do the forward transform (:ref:`fun-ispht`) 
   efficiently and I do that preprocessing from within :ref:`fun-scalarpatt`; second, putting the NumPy array within the :ref:`fun-scalarpatt`
   object makes the code consistent with how I deal with the :ref:`fun-vectorpatt` object, which has two NumPy arrays within it (:samp:`theta` and :samp:`phi`).
     
   


*spht* of My Head
=================


*spht* of Earth Image
=====================

.. note::
   See I did something here


