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
@newField description: Description
@newField revisions: Revisions
@newField departments: Departments
@newField applications: Applications

@author:
    rdireen

@organization:
    DireenTech Inc.

@description:

Credits:
The algorithms within this package have been implemented, in part, using the 
documentation within  the NIST Interagency/Internal Report (NISTIR) - 3955.
The majority or the code, however, has been developed using 
Ronald C. Wittmann's unpublished notes.

"""

#---------------------------------------------------------------------Built-ins
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import json
from os.path import dirname
from functools import wraps
import numbers

#---------------------------------------------------------------------3rd Party
import numpy as np

# TODO: Change all xrange instances to range
# and do a 'from six.moves import range' here
from six.moves import xrange  # @UnresolvedImport

#------------------------------------------------------------------------Custom

import spherepy.pysphi as pysphi  # python versions of the low level routines
import spherepy.ops as ops
import spherepy.csphi as csphi  # c extensions of the low level routines

#==============================================================================
# Global Declarations
#==============================================================================

with open(dirname(__file__) + '/pkg_info.json') as fp:
    _info = json.load(fp)
use_cext = _info["use_cext"]

scalar = 0
vector = 1

err_msg = {}
err_msg['nmax_g_mmax'] = "nmax must be greater than or equal to mmax"
err_msg['ncols_even'] = "the number of columns must be even"  
err_msg['tuple_2_el'] = "tuple must have 2 elements"
err_msg['ukn_coef_t'] = "unknown coefficient type"
err_msg['ukn_patt_t'] = "unknown pattern type"
err_msg['uknown_typ'] = "unknown type"
err_msg['m_out_bound'] = "m value out of bounds"
err_msg['n_out_bound'] = "n value out of bounds"
err_msg['no_v0_mode'] = "monopoles (n = 0) do not exist for type VectorCoefs"
err_msg['dim_no_mtch'] = "dimension mismatch"
err_msg['slice_noper'] = "slicing operation not permitted"
err_msg['idx_no_rec'] = "indexing method not recognized"
err_msg['SC_sz_msmtch'] = "ScalarCoefs could not be combined with sizes " + \
                          " (%d,%d) and (%d,%d)"
err_msg['VC_sz_msmtch'] = "VectorCoefs could not be combined with sizes " + \
                          "(%d,%d) and (%d,%d)"
err_msg['SP_sz_msmtch'] = "ScalarPatternUniform could not be combined " + \
                          " with sizes (%d,%d) and (%d,%d)"
err_msg['VP_sz_msmtch'] = "VectorPatternUniform could not be combined " + \
                          "with sizes (%d,%d) and (%d,%d)"
err_msg['no_combi_SC'] = "cannot combine type with ScalarCoefs"
err_msg['no_combi_VC'] = "cannot combine type with VectorCoefs"
err_msg['no_combi_SP'] = "cannot combine type with ScalarPatternUniform"
err_msg['no_combi_VP'] = "cannot combine type with VectorPatternUniform"
err_msg['vec1_s_vec2'] = "shape of vec1 and vec2 must be the the same"
err_msg['set_vc_val'] = "set coefficients like this vsc[n,m] = (val1,val2)"
err_msg['val_2tuple'] = "the value needs to be a tuple with length = 2 and" + \
                        " containing 2 numbers"
err_msg['set_sc1'] = "set scalar coefficients individually " + \
                     "(e.g. vcoef.scoef2[:,3] = vec)"
err_msg['set_sc2'] = "set scalar coefficients individually " + \
                                "(e.g. vcoef.scoef2[5,:] = vec)"
err_msg['shape_tc_pc'] = "shape of tcdata and pcdata must be the the same"
err_msg['nmax_too_lrg'] = "nmax value must be less than patterns nrows"
err_msg['mmax_too_lrg'] = "mmax value must be less than patterns ncols"
err_msg['use_mag_instead'] = "use mag( x ) instead"

err_msg['scoef_size'] = "vec must have length = " + \
                              "nmax + 1 + mmax * (2 * (nmax + 1) - mmax - 1)"

err_msg['vcoef_size'] = "vec1 and vec2 must have length = " + \
                              "nmax + 1 + mmax * (2 * (nmax + 1) - mmax - 1)"

err_msg['inverse_terr'] = "can't have nrows < nmax + 2 or " + \
                               "ncols < 2 * mmax + 2"

pretty_display_string = """

c[n, m]
=======

2: {7}  {4}  {2}  {6}  {8} 
1:                {3}  {1}  {5} 
0:                               {0}   
n  -------------  -------------  -------------  -------------  -------------  
       m = -2         m = -1         m = 0          m = 1          m = 2        
"""
                                
                                
#==============================================================================
# Classes
#==============================================================================

class ScalarCoefs(object):
    """This is the object for holding scalar spherical harmonic coefficients.
    To generate one of these objects it's best to used one of the factory 
    methods (e.g. spherepy.zeros_coefs, spherepy.random_coefs) or transform
    a pattern (ScalarPatternUniform) into a set of coefficients using 
    spherepy.spht. 

    Args:
          vec (numpy array complex128): Vector containing the coefficients.

          nmax (int): Largest *n* value in the set of modes.

          mmax (int): Largest abs(*m*) value in the set of modes.

    Example 1:
        Create a set of random coefficients::

            >>> c = spherepy.random_coefs(5,5)

    Example 2:
        Transform a pattern into coefficients::

            >>> f = spherepy.random_patt_uniform(6,8)
            >>> c = spherepy.spht(2,2)

    Example 2:
        Access *nmax* and *mmax*:

            >>> c.nmax
            >>> c.mmax

    **What are nmax and mmax?**

    Below is a diagram of a ScalarCoefs object with *nmax* = 5 and *mmax* = 2.
    If you count, you will see that there are 24 coefficients within *c*. 
    *nmax* is the largest *n* value in the object and *mmax* is the largest 
    *m* value in the structure (-*mmax* is the largest negative *m* value)::

        
        c[n, m]
        =======

        nmax = 5      *     *     *     *     *
        4:            *     *     *     *     *
        3:            *     *     *     *     *
        2:            *     *     *     *     *
        1:                  *     *     *
        0:                        *  
        n            ---   ---   ---   ---   ---
             -mmax = m=-2  m=-1  m=0   m=1   m=2 = mmax    
        

    .. note::
        In almost all cases, you can choose to set mmax to nmax. 

    """
    def __init__(self, vec, nmax, mmax):
        
        if mmax > nmax:
            raise ValueError(err_msg['nmax_g_mmax'])

        N = nmax + 1;
        NC = N + mmax * (2 * N - mmax - 1);
        if NC != len(vec):
            raise ValueError(err_msg['scoef_size'])

        self._vec = vec
        self._nmax = nmax
        self._mmax = mmax

    @property
    def nmax(self):
        """Largest *n* value."""
        return self._nmax

    @property
    def mmax(self):
        """Largest abs(*m*) value."""
        return self._mmax

    @property
    def size(self):
        """Total number of coefficients in the ScalarCoefs structure.

        Example::

            >>> sz  = c.size
            >>> N = c.nmax + 1
            >>> L = N+ c.mmax * (2 * N - c.mmax - 1);
            >>> assert sz == L
        """
        N = self.nmax + 1;
        NC = N + self.mmax * (2 * N - self.mmax - 1);
        assert NC == len(self._vec)
        return NC

    def copy(self):
        """Make a deep copy of this object.
        
        Example::

            >>> c2 = c.copy()
        
        """
        vec = np.copy(self._vec)
        return ScalarCoefs(vec, self.nmax, self.mmax)

    def window(self, vec):
        """Apply a window to the coefficients defined by *vec*. *vec* must
        have length *nmax* + 1.  This is good way to filter the pattern by
        windowing in the coefficient domain.

        Example::

            >>> vec = numpy.linspace(0, 1, c.nmax + 1)
            >>> c.window(vec)

        Args:
          vec (numpy.array): Vector of values to apply in the n direction of
          the data. Has length *nmax* + 1.

        Returns:
          Nothing, applies the window to the data in place.

        """

        slce = slice(None, None, None)
        
        self.__setitem__((slce, 0), self.__getitem__((slce, 0)) * vec)  
        for m in xrange(1, self.mmax + 1):
            self.__setitem__((slce, -m), self.__getitem__((slce, -m)) * vec[m:])
            self.__setitem__((slce, m), self.__getitem__((slce, m)) * vec[m:])

    def angular_power_spectrum(self):
        """Returns the angular power spectrum for the set of coefficients.
        That is, we compute

                   n
            c_n = sum  cnm * conj( cnm )
                  m=-n 

        Returns:
          power_spectrum (numpy.array, dtype=double) spectrum as a function of n.
        """

        # Added this routine as a result of my discussions with Ajinkya Nene	 
        #https://github.com/anene
        list_of_modes = self._reshape_m_vecs() 
        Nmodes = len(list_of_modes)

        angular_power = np.zeros( Nmodes, dtype = np.double)

        for n in range(0, Nmodes):
            mode = np.array( list_of_modes[n], dtype = np.complex128 )
            angular_power[n] = np.sum( np.abs(mode) ** 2 )
            
        return angular_power 

    def power(self):
        """ Returns the total power of the spectrum.
        That is, we compute

                Nmax  n
            P = sum  sum  cnm * conj( cnm )
                n=0  m=-n

        Returns:
           power (float): The total power in the coefficients.

        """

        # Added this routine as a result of my discussions with Ajinkya Nene	
        #https://github.com/anene
        return np.sum( np.abs(self._vec) ** 2)


    def _array_2d_repr(self):
        """creates a 2D array that has nmax + 1 rows and 2*mmax + 1 columns
        and provides a representation for the coefficients that makes 
        plotting easier"""

        sc_array = np.zeros((self.nmax + 1, 2 * self.mmax + 1),
                            dtype=np.complex128)

        lst = self._reshape_n_vecs()
        sc_array[0:self.nmax + 1, self.mmax] = lst[0]
        for m in xrange(1, self.mmax + 1):
            sc_array[m:self.nmax + 1, self.mmax - m] = lst[2 * m - 1]
            sc_array[m:self.nmax + 1, self.mmax + m] = lst[2 * m]

        return sc_array

    def _reshape_n_vecs(self):
        """return list of arrays, each array represents a different m mode"""

        lst = []
        sl = slice(None, None, None)
        lst.append(self.__getitem__((sl, 0)))
        for m in xrange(1, self.mmax + 1):
            lst.append(self.__getitem__((sl, -m)))
            lst.append(self.__getitem__((sl, m)))
        return lst

    def _reshape_m_vecs(self):
        """return list of arrays, each array represents a different n mode"""
        
        lst = []
        for n in xrange(0, self.nmax + 1):
            mlst = []
            if n <= self.mmax:
                nn = n
            else:
                nn = self.mmax            
            for m in xrange(-nn, nn + 1):
                mlst.append(self.__getitem__((n, m)))
            lst.append(mlst)
        return lst

    def __repr__(self):
        rs = self._reshape_m_vecs()
        return "ScalarCoefs((nmax = {0}, mmax = {1}) , {2} )".format(self.nmax,
                                                                     self.mmax,
                                                                     str(rs))

    def __str__(self):
        rs = self._reshape_m_vecs()
        st = ["ScalarCoefs(nmax = {0}, mmax = {1}):".format(self.nmax,
                                                          self.mmax)]
        for n, x in enumerate(rs):
            st.append("n = {0} :: {1}".format(n, x))
                      
        return '\n'.join(st)

    def __setitem__(self, arg, val):
        if len(arg) != 2:
            raise AttributeError(err_msg['tuple_2_el'])

        if isinstance(arg[0], slice) and isinstance(arg[1], int):

            m = arg[1]

            if ((arg[0].start is None) and 
                (arg[0].step is None) and 
                (arg[0].stop is None)):

                if np.abs(m) > self._mmax:
                    raise AttributeError(err_msg['m_out_bound'])

                idx_start = pysphi.mindx(m, self.nmax, self.mmax)
                idx_stop = idx_start + self.nmax - np.abs(m) + 1

                L = len(self._vec[idx_start:idx_stop])
                if len(val) != L:
                    raise AttributeError(err_msg['dim_no_mtch'])
                else:
                    self._vec[idx_start:idx_stop] = val
            else:
                raise AttributeError(err_msg['slice_noper'])

        elif isinstance(arg[0], int) and isinstance(arg[1], slice):
            
            n = arg[0]

            if ((arg[1].start is None) and 
                (arg[1].step is None) and 
                (arg[1].stop is None)):

                if (n > self.nmax) or (n < 0):
                    raise AttributeError(err_msg['n_out_bound'])

                if (len(val) != 2 * n + 1):
                    raise AttributeError(err_msg['dim_no_mtch'])

                if n <= self.mmax:
                    nn = n
                else:
                    nn = self.mmax            
                for k, m in enumerate(xrange(-nn, nn + 1)):
                    self.__setitem__((n, m), val[k])

            else:
                raise AttributeError(err_msg['slice_noper'])

        elif isinstance(arg[0], int) and isinstance(arg[1], int):

            n = arg[0]
            m = arg[1]

            if (n < 0) or (n > self._nmax):
                raise AttributeError(err_msg['n_out_bound'])

            if (np.abs(m) > n) or (np.abs(m) > self._mmax):
                raise AttributeError(err_msg['m_out_bound'])

            idx = pysphi.mindx(m, self.nmax, self.mmax) + n - np.abs(m)

            self._vec[idx] = val

        else:
            raise AttributeError(err_msg['idx_no_rec'])

    def __getitem__(self, arg):

        if len(arg) != 2:
            raise AttributeError(err_msg['tuple_2_el'])  

        if isinstance(arg[0], slice) and isinstance(arg[1], int):

            m = arg[1]

            if ((arg[0].start is None) and
               (arg[0].step is None) and
               (arg[0].stop is None)):

                if np.abs(m) > self._mmax:
                    raise AttributeError(err_msg['m_out_bound'])

                idx_start = pysphi.mindx(m, self.nmax, self.mmax)
                idx_stop = idx_start + self.nmax - np.abs(m) + 1

                return self._vec[idx_start:idx_stop]


            else:
                raise AttributeError(err_msg['slice_noper'])

        elif isinstance(arg[0], int) and isinstance(arg[1], slice):
            
            n = arg[0]

            if ((arg[1].start is None) and
               (arg[1].step is None) and
               (arg[1].stop is None)):

                if (n > self.nmax) or (n < 0):
                    raise AttributeError(err_msg['m_out_bound'])

                mlst = []
                if n <= self.mmax:
                    nn = n
                else:
                    nn = self.mmax            
                for m in xrange(-nn, nn + 1):
                    mlst.append(self.__getitem__((n, m)))

                return np.array(mlst, dtype=np.complex128)
            else:
                raise AttributeError(err_msg['slice_noper'])

        elif isinstance(arg[0], int) and isinstance(arg[1], int):

            n = arg[0]
            m = arg[1]

            if (n < 0) or (n > self._nmax):
                raise AttributeError(err_msg['n_out_bound'])

            if (np.abs(m) > n) or (np.abs(m) > self._mmax):
                raise AttributeError(err_msg['m_out_bound'])

            idx = pysphi.mindx(m, self.nmax, self.mmax) + n - np.abs(m)

            return self._vec[idx]

        elif isinstance(arg[0], slice) and isinstance(arg[1], slice):

            if ((arg[1].start is None) and
                (arg[1].step is None) and
                (arg[1].stop is None) and
                (arg[0].start == 0) and
                (arg[0].step is None) and
                isinstance(arg[0].stop, int)):
                
                if (arg[0].stop > self.nmax):
                    raise AttributeError(err_msg['n_out_bound'])

                new_nmax = arg[0].stop
                if new_nmax > self.mmax:
                    new_mmax = self.mmax
                else:
                    new_mmax = new_nmax

                N = new_nmax + 1;
                NC = N + new_mmax * (2 * N - new_mmax - 1);
                vec = np.zeros(NC, dtype=np.complex128)

                vec[0:new_nmax + 1] = self._vec[0:new_nmax + 1]
                for m in xrange(1, new_mmax + 1):
                    idx = pysphi.mindx(-m, self.nmax, self.mmax)
                    nidx = pysphi.mindx(-m, new_nmax, new_mmax)
                    ln = new_nmax - np.abs(m) + 1
                    vec[nidx:nidx + ln] = self._vec[idx:idx + ln]

                    idx = pysphi.mindx(m, self.nmax, self.mmax)
                    nidx = pysphi.mindx(m, new_nmax, new_mmax)
                    ln = new_nmax - np.abs(m) + 1
                    vec[nidx:nidx + ln] = self._vec[idx:idx + ln]

                return ScalarCoefs(vec, new_nmax, new_mmax)
            else:
                raise AttributeError(err_msg['idx_no_rec'])
        else:
            raise AttributeError(err_msg['idx_no_rec'])

    def _scalar_coef_op_left(func):
        """decorator for operator overloading when ScalarCoef is on the
        left"""
        @wraps(func)
        def verif(self, scoef):
            if isinstance(scoef, ScalarCoefs):
                if len(self._vec) == len(scoef._vec):
                    return ScalarCoefs(func(self, self._vec, scoef._vec),
                                        self.nmax,
                                        self.mmax)
                else:
                    raise ValueError(err_msg['SC_sz_msmtch'] % \
                                    (self.nmax, self.mmax,
                                     scoef.nmax, scoef.mmax))
        
            elif isinstance(scoef, numbers.Number):
                return ScalarCoefs(func(self, self._vec, scoef), self.nmax,
                                   self.mmax)
            else:
                raise TypeError(err_msg['no_combi_SC'])
        return verif

    def _scalar_coef_op_right(func):
        """decorator for operator overloading when ScalarCoef is on the
        right"""
        @wraps(func)
        def verif(self, scoef):
            if isinstance(scoef, numbers.Number):
                return ScalarCoefs(func(self, self._vec, scoef),
                                   self.nmax, self.mmax)
            else:
                raise TypeError(err_msg['no_combi_SC'])
        return verif

    @_scalar_coef_op_left
    def __add__(self, a, b):
        return a + b

    @_scalar_coef_op_right
    def __radd__(self, a, b):
        return  b + a

    @_scalar_coef_op_left
    def __sub__(self, a, b):
        return a - b

    @_scalar_coef_op_right
    def __rsub__(self, a, b):
        return b - a

    @_scalar_coef_op_left
    def __mul__(self, a, b):
        return a * b

    @_scalar_coef_op_right
    def __rmul__(self, a, b):
        return b * a

    @_scalar_coef_op_left
    def __div__(self, a, b):
        if isinstance(b, numbers.Number):
            if b == 0:
                return ZeroDivisionError()
            return a / b
        else:
            return a / b

    @_scalar_coef_op_right
    def __rdiv__(self, a, b):
        if isinstance(a, numbers.Number):
            if a == 0:
                return ZeroDivisionError()
            return b / a
        else:
            return b / a

    @_scalar_coef_op_left
    def __truediv__(self, a, b):
        if isinstance(b, numbers.Number):
            if b == 0:
                return ZeroDivisionError()
            return a / b
        else:
            return a / b

    @_scalar_coef_op_right
    def __rtruediv__(self, a, b):
        if isinstance(a, numbers.Number):
            if a == 0:
                return ZeroDivisionError()
            return b / a
        else:
            return b / a  


class VectorCoefs(object):
    """Holds the vector coefficients that represent a vector spherical pattern. 
    The function vspht returns this object"""
    def __init__(self, vec1, vec2, nmax, mmax):
        
        if(vec1.shape != vec2.shape):
            raise ValueError(err_msg['vec1_s_vec2'])
            
        if mmax > nmax:
            raise ValueError(err_msg['nmax_g_mmax'])

        N = nmax + 1;
        NC = N + mmax * (2 * N - mmax - 1);
        if NC != len(vec1):
            raise ValueError(err_msg['vcoef_size'])

        # There are no monopoles for structure VectorCoefs
        vec1[0] = 0 + 1j*0
        vec2[0] = 0 + 1j*0

        self.scoef1 = ScalarCoefs(vec1, nmax, mmax)
        self.scoef2 = ScalarCoefs(vec2, nmax, mmax)
        self._nmax = nmax
        self._mmax = mmax

    @property
    def nmax(self):
        """Largest n value."""
        return self._nmax

    @property
    def mmax(self):
        """Largest abs(m) value."""
        return self._mmax

    @property
    def size(self):
        """The number of modes in the structure"""
        N = self.nmax + 1;
        NC = N + self.mmax * (2 * N - self.mmax - 1);
        assert NC == len(self.scoef1._vec)
        assert NC == len(self.scoef2._vec)
        return NC

    def copy(self):
        """Make a deep copy of this object.
        
        Example::

            >>> c2 = c.copy()
        
        """
        vec1 = np.copy(self.scoef1._vec)
        vec2 = np.copy(self.scoef2._vec)
        return VectorCoefs(vec1, vec2, self.nmax, self.mmax)

    def _array_2d_repr(self):
        """creates a 2D array that has nmax + 1 rows and 2*mmax + 1 columns
        and provides a representation for the coefficients that makes 
        plotting easier"""

        return (self.scoef1._array_2d_repr(),
                self.scoef2._array_2d_repr())

    def _reshape_n_vecs(self):
        """return list of arrays, each array represents a different m mode"""

        return (self.scoef1._reshape_n_vecs(),
                self.scoef2._reshape_n_vecs())

    def _reshape_m_vecs(self):
        """return list of arrays, each array represents a different n mode"""
        
        return (self.scoef1._reshape_m_vecs(),
                self.scoef2._reshape_m_vecs())

    def __repr__(self):
        return "VectorCoefs(nmax = {0}, mmax = {1})".format(self.nmax,
                                                                     self.mmax)

    def __setitem__(self, arg, val):
        
        if len(arg) != 2:
            raise AttributeError(err_msg['tuple_2_el'])
        
        if isinstance(val, (list, tuple)):
            if len(val) != 2:
                raise AttributeError(err_msg['val_2tuple'])
        else:
            raise AttributeError(err_msg['val_2tuple'])

        if isinstance(arg[0], slice) and isinstance(arg[1], int):
            raise AttributeError(err_msg['set_sc1'])

        elif isinstance(arg[0], int) and isinstance(arg[1], slice):
            if arg[0] == 0:
                raise AttributeError(err_msg['no_v0_mode'])
            raise AttributeError(err_msg['set_sc2'])

        elif isinstance(arg[0], int) and isinstance(arg[1], int):
            n = arg[0]
            if n == 0:
                raise AttributeError(err_msg['no_v0_mode'])
            m = arg[1]

            try:
                self.scoef1[n, m] = val[0]
                self.scoef2[n, m] = val[1]
            except IndexError:
                AttributeError(err_msg['set_vc_val'])

        else:
            raise AttributeError(err_msg['idx_no_rec'])

    def __getitem__(self, arg):

        if len(arg) != 2:
            raise AttributeError(err_msg['tuple_2_el'])  

        if isinstance(arg[0], slice) and isinstance(arg[1], int):

            m = arg[1]

            if ((arg[0].start is None) and
               (arg[0].step is None) and
               (arg[0].stop is None)):

                return (self.scoef1[:, m], self.scoef2[:, m])
            else:
                raise AttributeError(err_msg['slice_noper'])

        elif isinstance(arg[0], int) and isinstance(arg[1], slice):
            
            n = arg[0]
            if n == 0:
                raise AttributeError(err_msg['no_v0_mode'])

            if ((arg[1].start is None) and
               (arg[1].step is None) and
               (arg[1].stop is None)):
                
                return (self.scoef1[n, :], self.scoef2[n, :])
            
            else:
                raise AttributeError(err_msg['slice_noper'])

        elif isinstance(arg[0], int) and isinstance(arg[1], int):
            n = arg[0]
            if n == 0:
                raise AttributeError(err_msg['no_v0_mode'])
            m = arg[1]
            return (self.scoef1[n, m], self.scoef2[n, m])

        elif isinstance(arg[0], slice) and isinstance(arg[1], slice):

            if ((arg[1].start is None) and
                (arg[1].step is None) and
                (arg[1].stop is None) and
                (arg[0].start == 0) and
                (arg[0].step is None) and
                isinstance(arg[0].stop, int)):
                
                if (arg[0].stop > self.nmax):
                    raise AttributeError(err_msg['n_out_bound'])

                new_nmax = arg[0].stop
                if new_nmax > self.mmax:
                    new_mmax = self.mmax
                else:
                    new_mmax = new_nmax
                return VectorCoefs(self.scoef1[0:new_nmax, :]._vec,
                                   self.scoef2[0:new_nmax, :]._vec,
                                   new_nmax, new_mmax)
            else:
                raise AttributeError(err_msg['idx_no_rec'])

        else:
            raise AttributeError(err_msg['idx_no_rec'])

    def _vector_coef_op_left(func):
        """decorator for operator overloading when VectorCoef is on the
        left"""
        @wraps(func)
        def verif(self, vcoef):
            if isinstance(vcoef, VectorCoefs):
                if len(vcoef.scoef1._vec) == len(vcoef.scoef1._vec):
                    return VectorCoefs(func(self, self.scoef1._vec,
                                                  vcoef.scoef1._vec),
                                       func(self, self.scoef2._vec,
                                                  vcoef.scoef2._vec),
                                        self.nmax,
                                        self.mmax)
                else:
                    raise ValueError(err_msg['VC_sz_msmtch'] % \
                                    (self.nmax, self.mmax,
                                     vcoef.nmax, vcoef.mmax))
        
            elif isinstance(vcoef, numbers.Number):
                return VectorCoefs(func(self, self.scoef1._vec, vcoef),
                                   func(self, self.scoef2._vec, vcoef),
                                   self.nmax,
                                   self.mmax)
            else:
                raise TypeError(err_msg['no_combi_VC'])
        return verif

    def _vector_coef_op_right(func):
        """decorator for operator overloading when VectorCoefs is on the
        right"""
        @wraps(func)
        def verif(self, vcoef):
            if isinstance(vcoef, numbers.Number):
                return VectorCoefs(func(self, self.scoef1._vec, vcoef),
                                   func(self, self.scoef2._vec, vcoef),
                                   self.nmax, self.mmax)
            else:
                raise TypeError(err_msg['no_combi_VC'])
        return verif

    @_vector_coef_op_left
    def __add__(self, a, b):
        return a + b

    @_vector_coef_op_right
    def __radd__(self, a, b):
        return  b + a

    @_vector_coef_op_left
    def __sub__(self, a, b):
        return a - b

    @_vector_coef_op_right
    def __rsub__(self, a, b):
        return b - a

    @_vector_coef_op_left
    def __mul__(self, a, b):
        return a * b

    @_vector_coef_op_right
    def __rmul__(self, a, b):
        return b * a

    @_vector_coef_op_left
    def __div__(self, a, b):
        if isinstance(b, numbers.Number):
            if b == 0:
                return ZeroDivisionError()
            return a / b
        else:
            return a / b

    @_vector_coef_op_right
    def __rdiv__(self, a, b):
        if isinstance(a, numbers.Number):
            if a == 0:
                return ZeroDivisionError()
            return b / a
        else:
            return b / a  

    @_vector_coef_op_left
    def __truediv__(self, a, b):
        if isinstance(b, numbers.Number):
            if b == 0:
                return ZeroDivisionError()
            return a / b
        else:
            return a / b

    @_vector_coef_op_right
    def __rtruediv__(self, a, b):
        if isinstance(a, numbers.Number):
            if a == 0:
                return ZeroDivisionError()
            return b / a
        else:
            return b / a  


class ScalarPatternUniform(object):
    """This is the object for holding a scalar pattern. 
    
    .. note::

        The scalar pattern must satisfy certain symmetry constraints for the
        *spht* function to operate on it. This is really the only reason why
        you wouldn't just pass a 2D NumPy array directly to *spht*. 

    Args:
          cdata (2D numpy array complex128): Pattern to be transformed.

          doublesphere (bool, optional): Don't worry about it.
    
    Example 1
        Creating a ScalarPatternUniform object from a NumPy array::

            >>> a = numpy.zeros((5, 8), dtype = numpy.complex128)
            >>> f = ScalarPatternUniform(a)

        The above will put the array into the proper form for processing 
        with the *spht*.

    Example 2:
        Same as Example 1, but using a factory I made::

            >>> f = spherepy.zeros_patt_uniform(5, 8)


    **How do I put my own pattern in?**

    You must have a full set of data over the surface of a sphere. The 
    theta coordinate runs *0* to *pi* radians from the north to south poles,
    and corresponds to the rows of the NumPy array. The phi coordinate runs
    from 0 to *2*pi - 2*pi / ncols*  radians, and corresponds to the colums
    of the NumPy array.
    
    As an example, consider a pattern *f* with *nrows* = 5 and *ncols* = 8. 
    If you place an array with these dimensions into ScalarPatternUniform, the
    code will assume each row and column corresponds to a position on the 
    sphere, as depicted in the following table::

        f[n_row, n_col]
        ===============

               
              : phi  0     pi/4   pi/2  3*pi/4  pi  5*pi/4  3*pi/2  7*pi/4
         theta      ---    ---    ---    ---    ---    ---    ---    ---
             0:      *      *      *      *      *      *      *      *
          pi/4:      *      *      *      *      *      *      *      *
          pi/2:      *      *      *      *      *      *      *      *
        3*pi/4:      *     5.0     *      *      *      *      *      *
            pi:      *      *      *      *      *      *      *      *  
       
    Therefore  *f[3, 1]*, which indexes the  value 5.0 in the array above, 
    corresponds to location *theta = 3*pi/4*, *phi = pi/4* on the sphere. 

    **Some notes:**

    It should be clear that *f[0, :]* corresponds to the north pole of the 
    sphere, and *f[4, :]* corresponds to the south pole of the sphere.

    .. note::
        **BIG NOTE:** I think the most difficult part of using this code is 
        getting the pattern data in right. If you're getting strange 
        results, come back to this section often to make sure you understand
        how the data is supposed to be entered. 

        
   

    """

    def __init__(self, cdata, doublesphere=False):

        if(doublesphere == False):
            self._dsphere = continue_sphere(cdata, 1)
        else:
            self._dsphere = cdata

    def __repr__(self):

        return "ScalarPatternUniform({0})".format(self.array.__repr__())

    @property
    def nrows(self):
        """Returns the number of rows in the NumPy array."""

        return int(self._dsphere.shape[0] / 2) + 1

    @property
    def ncols(self):
        """Returns the number of columns in the NumPy array."""

        return self._dsphere.shape[1]

    @property
    def shape(self):
        """Same as patt.cdata.shape"""

        return (self.nrows, self.ncols)

    @property
    def doublesphere(self):
        """ Returns a NumPy array with the proper symmetry for 
        transformation."""

        return self._dsphere

    @property
    def array(self):
        """Returns the NumPy array (same as cdata)"""

        return self._dsphere[0:self.nrows, :]

    @property
    def cdata(self):
        """Return the NumPy array."""

        return self._dsphere[0:self.nrows, :]

    @property
    def is_symmetric(self):
        """return true if the data has the proper symmetry for transformation.
        """
        if self.single_val < 1e-14:
            return True 
        else:
            return False

    @property
    def single_val(self):
        """return relative error of worst point that might make the data none
        symmetric.
        """
        nrows = self._dsphere.shape[0]
        ncols = self._dsphere.shape[1]

        sv = 0.0
        for n in xrange(0, nrows):
            for m in xrange(0, ncols):
                s = self._dsphere[np.mod(nrows - n, nrows),
                              np.mod(int(np.floor(ncols / 2)) + m, ncols)]
                t = self._dsphere[n, m]

                if s != 0:
                    sabs = np.abs((s - t) / s)
                    if sabs >= sv:
                        sv = sabs
                elif t != 0:
                    sabs = 1.0
                    if sabs >= sv:
                        sv = sabs

        return sv

    def _scalar_pattern_uniform_op_left(func):
        """Decorator for operator overloading when ScalarPatternUniform is on 
        the left."""
        @wraps(func)
        def verif(self, patt):
            if isinstance(patt, ScalarPatternUniform):
                if self._dsphere.shape == patt._dsphere.shape:
                    return ScalarPatternUniform(func(self, self._dsphere,
                                                     patt._dsphere),
                                                doublesphere=True)
                else:
                    raise ValueError(err_msg['SP_sz_msmtch'] % \
                                            (self.nrows, self.ncols,
                                            patt.nrows, patt.ncols))
        
            elif isinstance(patt, numbers.Number):
                return ScalarPatternUniform(func(self, self._dsphere, patt),
                                            doublesphere=True)
            else:
                raise TypeError(err_msg['no_combi_SP'])
        return verif

    def _scalar_pattern_uniform_op_right(func):
        """Decorator for operator overloading when ScalarPatternUniform is on
        the right."""
        @wraps(func)
        def verif(self, patt):
            if isinstance(patt, numbers.Number):
                return ScalarPatternUniform(func(self, self._dsphere, patt),
                                   doublesphere=True)
            else:
                raise TypeError(err_msg['no_combi_SP'])
        return verif

    @_scalar_pattern_uniform_op_left
    def __add__(self, a, b):
        return a + b

    @_scalar_pattern_uniform_op_right
    def __radd__(self, a, b):
        return  b + a

    @_scalar_pattern_uniform_op_left
    def __sub__(self, a, b):
        return a - b

    @_scalar_pattern_uniform_op_right
    def __rsub__(self, a, b):
        return b - a

    @_scalar_pattern_uniform_op_left
    def __mul__(self, a, b):
        return a * b

    @_scalar_pattern_uniform_op_right
    def __rmul__(self, a, b):
        return b * a

    @_scalar_pattern_uniform_op_left
    def __div__(self, a, b):
        if isinstance(b, numbers.Number):
            if b == 0:
                return ZeroDivisionError()
            return a / b
        else:
            return a / b

    @_scalar_pattern_uniform_op_right
    def __rdiv__(self, a, b):
        if isinstance(a, numbers.Number):
            if a == 0:
                return ZeroDivisionError()
            return b / a
        else:
            return b / a 

    @_scalar_pattern_uniform_op_left
    def __truediv__(self, a, b):
        if isinstance(b, numbers.Number):
            if b == 0:
                return ZeroDivisionError()
            return a / b
        else:
            return a / b

    @_scalar_pattern_uniform_op_right
    def __rtruediv__(self, a, b):
        if isinstance(a, numbers.Number):
            if a == 0:
                return ZeroDivisionError()
            return b / a
        else:
            return b / a  


class ScalarPatternNonUniform:
    pass

class TransversePatternUniform:
    
    def __init__(self, tcdata, pcdata, doublesphere=False):
        
        if(tcdata.shape != pcdata.shape):
            raise ValueError(err_msg['shape_tc_pc'])

        if(doublesphere == False):
            self._tdsphere = continue_sphere(tcdata, -1)
            self._pdsphere = continue_sphere(pcdata, -1)
        else:
            self._tdsphere = tcdata
            self._pdsphere = pcdata

    def __repr__(self):

        return self.array.__repr__()

    @property
    def nrows(self):
        return int(self._tdsphere.shape[0] / 2) + 1 

    @property
    def ncols(self):
        return self._tdsphere.shape[1]

    @property
    def shape(self):
        return (self.nrows, self.ncols)

    @property
    def doublesphere(self):
        return (self._tdsphere, self._pdsphere)

    @property
    def array(self):
        return (self._tdsphere[0:self.nrows, :],
                self._pdsphere[0:self.nrows, :])

    @property
    def theta(self):
        return self._tdsphere[0:self.nrows, :]

    @property
    def phi(self):
        return self._pdsphere[0:self.nrows, :]

    @property
    def theta_double(self):
        return self._tdsphere

    @property
    def phi_double(self):
        return self._pdsphere

    @property
    def is_symmetric(self):
        """return true if the data is symmetric"""
        if (self.single_val[0] < 1e-13) and (self.single_val[1] < 1e-13):
            return True 
        else:
            return False

    @property
    def single_val(self):
        """return relative error of worst point that might make the data none
        symmetric.
        """
        
        sv_t = self._sv(self._tdsphere)
        sv_p = self._sv(self._tdsphere)
        
        return (sv_t, sv_p)
        
    
    def _sv(self, dat):
        
        nrows = dat.shape[0]
        ncols = dat.shape[1]

        sv = 0.0
        for n in xrange(0, nrows):
            for m in xrange(0, ncols):
                s = -dat[np.mod(nrows - n, nrows),
                              np.mod(int(np.floor(ncols / 2)) + m, ncols)]
                t = dat[n, m]

                if s != 0:
                    sabs = np.abs((s - t) / s)
                    if sabs >= sv:
                        sv = sabs
                elif t != 0:
                    sabs = 1.0
                    if sabs >= sv:
                        sv = sabs

        return sv
        
    def _vector_pattern_uniform_op_left(func):
        """decorator for operator overloading when VectorPatternUniform is on 
        the left"""
        @wraps(func)
        def verif(self, patt):
            if isinstance(patt, TransversePatternUniform):
                if self._tdsphere.shape == patt._tdsphere.shape:
                    return TransversePatternUniform(func(self, self._tdsphere,
                                                     patt._tdsphere),
                                                func(self, self._pdsphere,
                                                     patt._pdsphere),
                                                doublesphere=True)
                else:
                    raise ValueError(err_msg['VP_sz_msmtch'] % \
                                            (self.nrows, self.ncols,
                                            patt.nrows, patt.ncols))
        
            elif isinstance(patt, numbers.Number):
                return TransversePatternUniform(func(self, self._tdsphere, patt),
                                            func(self, self._pdsphere, patt),
                                            doublesphere=True)
            else:
                raise TypeError(err_msg['no_combi_VP'])
        return verif

    def _vector_pattern_uniform_op_right(func):
        """decorator for operator overloading when VectorPatternUniform is on
        the right"""
        @wraps(func)
        def verif(self, patt):
            if isinstance(patt, numbers.Number):
                return TransversePatternUniform(func(self, self._tdsphere, patt),
                                            func(self, self._pdsphere, patt),
                                            doublesphere=True)
            else:
                raise TypeError(err_msg['no_combi_VP'])
        return verif

    @_vector_pattern_uniform_op_left
    def __add__(self, a, b):
        return a + b

    @_vector_pattern_uniform_op_right
    def __radd__(self, a, b):
        return  b + a

    @_vector_pattern_uniform_op_left
    def __sub__(self, a, b):
        return a - b

    @_vector_pattern_uniform_op_right
    def __rsub__(self, a, b):
        return b - a

    @_vector_pattern_uniform_op_left
    def __mul__(self, a, b):
        return a * b

    @_vector_pattern_uniform_op_right
    def __rmul__(self, a, b):
        return b * a

    @_vector_pattern_uniform_op_left
    def __div__(self, a, b):
        if isinstance(b, numbers.Number):
            if b == 0:
                return ZeroDivisionError()
            return a / b
        else:
            return a / b  

    @_vector_pattern_uniform_op_right
    def __rdiv__(self, a, b):
        if isinstance(a, numbers.Number):
            if a == 0:
                return ZeroDivisionError()
            return b / a
        else:
            return b / a 

    @_vector_pattern_uniform_op_left
    def __truediv__(self, a, b):
        if isinstance(b, numbers.Number):
            if b == 0:
                return ZeroDivisionError()
            return a / b
        else:
            return a / b  

    @_vector_pattern_uniform_op_right
    def __rtruediv__(self, a, b):
        if isinstance(a, numbers.Number):
            if a == 0:
                return ZeroDivisionError()
            return b / a
        else:
            return b / a    


class TransversePatternNonUniform:
    pass

#==============================================================================
# Functions
#==============================================================================
        
#---------------------------------------------------------------------Factories
def zeros_coefs(nmax, mmax, coef_type=scalar):
    """Returns a ScalarCoefs object or a VectorCoeffs object where each of the 
    coefficients is set to 0. The structure is such that *nmax* is th largest 
    *n* can be in c[n, m], and *mmax* is the largest *m* can be for any *n*.
    (See *ScalarCoefs* and *VectorCoefs* for details.)
    

    Examples::
        
        >>> c = spherepy.zeros_coefs(5, 3, coef_type = spherepy.scalar)
        >>> c = spherepy.zeros_coefs(5, 3) # same as above
        >>> vc = spherepy.zeros_coefs(5, 3, coef_type = spherepy.vector)

    Args:
          nmax (int): Largest *n* value in the set of modes.

          mmax (int): Largest abs(*m*) value in the set of modes.

          coef_type (int, optional): Set to 0 for scalar, and 1 for vector.
          The default option is scalar. If you would like to return a set of 
          vector spherical hamonic coefficients, the preferred way to do so 
          is vc = spherepy.zeros_coefs( 10, 12, coef_type = spherepy.vector).

    Returns:
      coefs: Returns a ScalarCoefs object if coef_type is either blank or
      set to 0. Returns a VectorCoefs object if coef_type = 1.

    Raises:
      TypeError: If coef_type is anything but 0 or 1.

    """

    if(mmax > nmax):
        raise ValueError(err_msg['nmax_g_mmax'])

    if(coef_type == scalar):
        L = (nmax + 1) + mmax * (2 * nmax - mmax + 1)
        vec = np.zeros(L, dtype=np.complex128)
        return  ScalarCoefs(vec, nmax, mmax)
    elif(coef_type == vector):
        L = (nmax + 1) + mmax * (2 * nmax - mmax + 1)
        vec1 = np.zeros(L, dtype=np.complex128)
        vec2 = np.zeros(L, dtype=np.complex128)
        return  VectorCoefs(vec1, vec2, nmax, mmax)
    else:
        raise TypeError(err_msg['ukn_coef_t'])

def ones_coefs(nmax, mmax, coef_type=scalar):
    """Returns a ScalarCoefs object or a VectorCoeffs object where each of the 
    coefficients is set to 1.0. The structure is such that *nmax* is th largest 
    *n* can be in c[n, m], and *mmax* is the largest *m* can be for any *n*.
    (See *ScalarCoefs* and *VectorCoefs* for details.)
    

    Examples::
        
        >>> c = spherepy.ones_coefs(5, 3, coef_type = spherepy.scalar)
        >>> c = spherepy.ones_coefs(5, 3) # same as above
        >>> vc = spherepy.ones_coefs(5, 3, coef_type = spherepy.vector)

    Args:
          nmax (int): Largest *n* value in the set of modes.

          mmax (int): Largest abs(*m*) value in the set of modes.

          coef_type (int, optional): Set to 0 for scalar, and 1 for vector.
          The default option is scalar. If you would like to return a set of 
          vector spherical hamonic coefficients, the preferred way to do so 
          is vc = spherepy.zeros_coefs( 10, 12, coef_type = spherepy.vector).

    Returns:
      coefs: Returns a ScalarCoefs object if coef_type is either blank or
      set to 0. Returns a VectorCoefs object if coef_type = 1.

    Raises:
      TypeError: If coef_type is anything but 0 or 1.

    """

    if(mmax > nmax):
        raise ValueError(err_msg['nmax_g_mmax'])
    
    if(coef_type == scalar):
        L = (nmax + 1) + mmax * (2 * nmax - mmax + 1)
        vec = np.ones(L, dtype=np.complex128)
        return  ScalarCoefs(vec, nmax, mmax)
    elif(coef_type == vector):
        L = (nmax + 1) + mmax * (2 * nmax - mmax + 1)
        vec1 = np.ones(L, dtype=np.complex128)
        vec1[0] = 0
        vec2 = np.ones(L, dtype=np.complex128)
        vec2[0] = 0
        return  VectorCoefs(vec1, vec2, nmax, mmax)
    else:
        raise TypeError(err_msg['ukn_coef_t'])

def random_coefs(nmax, mmax, mu=0.0, sigma=1.0, coef_type=scalar):
    """Returns a ScalarCoefs object or a VectorCoeffs object where each of the 
    coefficients is a normal random variable with mean 0 and standardard 
    deviation 1.0. The structure is such that *nmax* is th largest 
    *n* can be in c[n, m], and *mmax* is the largest *m* can be for any *n*.
    (See *ScalarCoefs* and *VectorCoefs* for details.)
    

    Examples::
        
        >>> c = spherepy.random_coefs(5, 3, coef_type = spherepy.scalar)
        >>> c = spherepy.random_coefs(5, 3) # same as above
        >>> vc = spherepy.random_coefs(5, 3, coef_type = spherepy.vector)

    Args:
          nmax (int): Largest *n* value in the set of modes.

          mmax (int): Largest abs(*m*) value in the set of modes.

          coef_type (int, optional): Set to 0 for scalar, and 1 for vector.
          The default option is scalar. If you would like to return a set of 
          vector spherical hamonic coefficients, the preferred way to do so 
          is vc = spherepy.zeros_coefs( 10, 12, coef_type = spherepy.vector).

    Returns:
      coefs: Returns a ScalarCoefs object if coef_type is either blank or
      set to 0. Returns a VectorCoefs object if coef_type = 1.

    Raises:
      TypeError: If coef_type is anything but 0 or 1.

    """

    if(mmax > nmax):
        raise ValueError(err_msg['nmax_g_mmax'])

    if(coef_type == scalar):
        L = (nmax + 1) + mmax * (2 * nmax - mmax + 1)
        vec = np.random.normal(mu, sigma, L) + \
              1j * np.random.normal(mu, sigma, L)
        return  ScalarCoefs(vec, nmax, mmax)
    elif(coef_type == vector):
        L = (nmax + 1) + mmax * (2 * nmax - mmax + 1)
        vec1 = np.random.normal(mu, sigma, L) + \
              1j * np.random.normal(mu, sigma, L)
        vec1[0] = 0
        vec2 = np.random.normal(mu, sigma, L) + \
              1j * np.random.normal(mu, sigma, L)
        vec2[0] = 0
        return  VectorCoefs(vec1, vec2, nmax, mmax)
    else:
        raise TypeError(err_msg['ukn_coef_t'])

def zeros_patt_uniform(nrows, ncols, patt_type=scalar):
    """Returns a ScalarPatternUniform object or a VectorPatternUniform object
    where each of the elements is set to 0. *nrows* is the number of rows in 
    the pattern, which corresponds to the theta axis. *ncols* must be even
    and is the number of columns in the pattern and corresponds to the phi 
    axis.
    (See *ScalarPatternUniform* and *VectorPatternUniform* for details.)
    
    Examples::
        
        >>> f = spherepy.zeros_patt_uniform(6, 8, coef_type = spherepy.scalar)
        >>> f = spherepy.zeros_patt_uniform(6, 8) # same as above
        >>> F = spherepy.zeros_patt_uniform(6, 8, coef_type = spherepy.vector)

    Args:
          nrows (int): Number of rows corresponding to the theta axis.

          ncols (int): Number of columns corresponding to the phi axis. To get 
          the speed and accuracy I need, this value **must** be even. 

          coef_type (int, optional): Set to 0 for scalar, and 1 for vector.
          The default option is scalar. 

    Returns:
      coefs: Returns a ScalarPatternUniform object if coef_type is either 
      blank or set to 0. Returns a VectorPatternUniform object if 
      coef_type = 1.

    Raises:
      ValueError: If ncols is not even.

      TypeError: If coef_type is anything but 0 or 1.

    """
    
    if np.mod(ncols, 2) == 1:
        raise ValueError(err_msg['ncols_even'])

    if(patt_type == scalar):
        cdata = np.zeros((2 * nrows - 2, ncols), dtype=np.complex128)
        return ScalarPatternUniform(cdata, doublesphere=True)

    elif(patt_type == vector):
        tcdata = np.zeros((2 * nrows - 2, ncols), dtype=np.complex128)
        pcdata = np.zeros((2 * nrows - 2, ncols), dtype=np.complex128)
        return TransversePatternUniform(tcdata, pcdata, doublesphere=True)

    else:
        raise TypeError(err_msg['ukn_patt_t'])

def ones_patt_uniform(nrows, ncols, patt_type=scalar):
    """Returns a ScalarPatternUniform object or a VectorPatternUniform object
    where each of the elements is set to 1. *nrows* is the number of rows in 
    the pattern, which corresponds to the theta axis. *ncols* must be even
    and is the number of columns in the pattern and corresponds to the phi 
    axis.
    (See *ScalarPatternUniform* and *VectorPatternUniform* for details.)
    
    Examples::
        
        >>> f = spherepy.ones_patt_uniform(6, 8, coef_type = spherepy.scalar)
        >>> f = spherepy.ones_patt_uniform(6, 8) # same as above
        >>> F = spherepy.ones_patt_uniform(6, 8, coef_type = spherepy.vector)

    Args:
          nrows (int): Number of rows corresponding to the theta axis.

          ncols (int): Number of columns corresponding to the phi axis. To get 
          the speed and accuracy I need, this value **must** be even. 

          coef_type (int, optional): Set to 0 for scalar, and 1 for vector.
          The default option is scalar. 

    Returns:
      coefs: Returns a ScalarPatternUniform object if coef_type is either 
      blank or set to 0. Returns a VectorPatternUniform object if 
      coef_type = 1.

    Raises:
      ValueError: If ncols is not even.

      TypeError: If coef_type is anything but 0 or 1.

    """
    
    if np.mod(ncols, 2) == 1:
        raise ValueError(err_msg['ncols_even'])

    if(patt_type == scalar):
        cdata = np.ones((2 * nrows - 2, ncols), dtype=np.complex128)
        return ScalarPatternUniform(cdata, doublesphere=True)

    elif(patt_type == vector):
        tcdata = np.ones((2 * nrows - 2, ncols), dtype=np.complex128)
        pcdata = np.ones((2 * nrows - 2, ncols), dtype=np.complex128)
        return TransversePatternUniform(tcdata, pcdata, doublesphere=True)

    else:
        raise TypeError(err_msg['ukn_patt_t'])

def random_patt_uniform(nrows, ncols, patt_type=scalar):
    """Returns a ScalarPatternUniform object or a VectorPatternUniform object
    where each of the elements is set to a normal random variable with zero 
    mean and unit standard deviation. *nrows* is the number of rows in 
    the pattern, which corresponds to the theta axis. *ncols* must be even
    and is the number of columns in the pattern and corresponds to the phi 
    axis.
    (See *ScalarPatternUniform* and *VectorPatternUniform* for details.)
    
    Examples::
        
        >>> f = spherepy.random_patt_uniform(6, 8, coef_type = spherepy.scalar)
        >>> f = spherepy.random_patt_uniform(6, 8) # same as above
        >>> F = spherepy.random_patt_uniform(6, 8, coef_type = spherepy.vector)

    Args:
          nrows (int): Number of rows corresponding to the theta axis.

          ncols (int): Number of columns corresponding to the phi axis. To get 
          the speed and accuracy I need, this value **must** be even. 

          coef_type (int, optional): Set to 0 for scalar, and 1 for vector.
          The default option is scalar. 

    Returns:
      coefs: Returns a ScalarPatternUniform object if coef_type is either 
      blank or set to 0. Returns a VectorPatternUniform object if 
      coef_type = 1.

    Raises:
      ValueError: If ncols is not even.

      TypeError: If coef_type is anything but 0 or 1.

    """

    if np.mod(ncols, 2) == 1:
        raise ValueError(err_msg['ncols_even'])

    if(patt_type == scalar):
        vec = np.random.normal(0.0, 1.0, nrows * ncols) + \
              1j * np.random.normal(0.0, 1.0, nrows * ncols)
        return ScalarPatternUniform(vec.reshape((nrows, ncols)),
                                    doublesphere=False)

    elif(patt_type == vector):
        vec1 = np.random.normal(0.0, 1.0, nrows * ncols) + \
              1j * np.random.normal(0.0, 1.0, nrows * ncols)
        vec2 = np.random.normal(0.0, 1.0, nrows * ncols) + \
              1j * np.random.normal(0.0, 1.0, nrows * ncols)
        return TransversePatternUniform(vec1.reshape((nrows, ncols)),
                                    vec2.reshape((nrows, ncols)),
                                    doublesphere=False)

    else:
        raise TypeError(err_msg['ukn_patt_t'])

def zeros_patt_nonuniform(theta_phi, patt_type=scalar):
    raise NotImplementedError()

def ones_patt_nonuniform(theta_phi, patt_type=scalar):
    raise NotImplementedError()

def random_patt_nonuniform(theta_phi, ncols, patt_type=scalar):
    raise NotImplementedError()

#------------------------------------------------------------------------------ 
def array(patt):

    if isinstance(patt, ScalarPatternUniform):
        return patt.array
    elif isinstance(patt, TransversePatternUniform):
        return patt.array
    else:
        raise TypeError(err_msg['uknown_typ'])
    
def abs(sobj):
    if isinstance(sobj, ScalarPatternUniform):
        return np.abs(sobj.array)
    elif isinstance(sobj, ScalarCoefs):
        return ScalarCoefs(np.abs(sobj._vec), sobj.nmax, sobj.mmax)
    else:
        raise TypeError(err_msg['uknown_typ'])

def mag2(sobj):
    if isinstance(sobj, ScalarPatternUniform):
        return np.abs(sobj.array) ** 2
    elif isinstance(sobj, ScalarCoefs):
        return ScalarCoefs(np.abs(sobj._vec), sobj.nmax, sobj.mmax)
    elif isinstance(sobj, TransversePatternUniform):
        return np.abs(sobj.theta) ** 2 + np.abs(sobj.phi) ** 2
    elif isinstance(sobj, VectorCoefs):
        return ScalarCoefs(np.abs(sobj.scoef1._vec) ** 2 + \
                           np.abs(sobj.scoef2._vec) ** 2,
                           sobj.nmax, sobj.mmax)
    else:
        raise TypeError(err_msg['uknown_typ'])

def mag(sobj):
    if isinstance(sobj, ScalarPatternUniform):
        return np.abs(sobj.array)
    elif isinstance(sobj, ScalarCoefs):
        return ScalarCoefs(np.abs(sobj._vec), sobj.nmax, sobj.mmax)
    elif isinstance(sobj, TransversePatternUniform):
        return np.sqrt(mag2(sobj))
    elif isinstance(sobj, VectorCoefs):
        return ScalarCoefs(np.sqrt(mag2(sobj)._vec),
                           sobj.nmax, sobj.mmax)
    else:
        raise TypeError(err_msg['uknown_typ'])
    
def compare_relative(sc1, sc2):  
    num = np.sum(np.abs(sc1._vec - sc2._vec))
    den = np.sum(np.abs(sc1._vec))
    return num / den
    
def continue_sphere(cdata, sym):
    
    nrows = cdata.shape[0]
    ncols = cdata.shape[1]

    zdata = np.zeros([nrows - 2, ncols], dtype=np.complex128)
    ddata = np.concatenate([cdata, zdata])
    return double_sphere(ddata, sym)


def double_sphere(cdata, sym):
    """ Ensures that the data within cdata has double sphere symmetry.

    Example::

        >>> spherepy.doublesphere(cdata, 1)

    Args:
        sym (int): is 1 for scalar data and -1 for vector data

    Returns:
        numpy.array([*,*], dtype=np.complex128) containing array with 
        doublesphere symmetry.
    """
    
    nrows = cdata.shape[0]
    ncols = cdata.shape[1]

    ddata = np.zeros([nrows, ncols], dtype=np.complex128)

    for n in xrange(0, nrows):
        for m in xrange(0, ncols):
            s = sym * cdata[np.mod(nrows - n, nrows),
                          np.mod(int(np.floor(ncols / 2)) + m, ncols)]
            t = cdata[n, m]

            if s * t == 0:
                ddata[n, m] = s + t
            else:
                ddata[n, m] = (s + t) / 2

    return ddata

def _tiny_rep(c):
    sr = "{0:.2}".format(c)
    if sr[0] == '(':
        sr = sr[1:-1]
    return sr

def pretty_coefs(c):
    """Prints out the first 2 modes of a ScalarCoeffs object. This is mostly 
    used for instructional purposes.
    (*ScalarPatternUniform*)
    
    Example::
        
        >>> spherepy.pretty_coefs(c)
        
        c[n, m]
        =======

        2:       0j             0j            0j             0j             0j 
        1:                     1+0j           0j           -1+0j 
        0:                                    1j  
        n  -------------  -------------  -------------  -------------  -------------  
               m = -2         m = -1         m = 0          m = 1          m = 2    

    Args:
          c (ScalarCoefs): Coefficients to be printed.

    Returns:
      nothing, just outputs something pretty to the console.

    """

    cfit = c[0:2, :]
    cvec = cfit._vec

    sa = [_tiny_rep(val) for val in cvec]

    while len(sa) < 9:
        sa.append("")

    sa = [sa[n].center(13) for n in range(0, 9)]

    print(pretty_display_string.format(sa[0], sa[1], sa[2],
                                       sa[3], sa[4], sa[5],
                                       sa[6], sa[7], sa[8]))


def spht(ssphere, nmax=None, mmax=None):
    """Transforms ScalarPatternUniform object *ssphere* into a set of scalar
    spherical harmonics stored in ScalarCoefs.

    Example::

        >>> p = spherepy.random_patt_uniform(6, 8)
        >>> c = spherepy.spht(p)
        >>> spherepy.pretty_coefs(c)

    Args:
      ssphere (ScalarPatternUniform): The pattern to be transformed.

      nmax (int, optional): The maximum number of *n* values required. If a 
      value isn't passed, *nmax* is the number of rows in ssphere minus one.

      mmax (int, optional): The maximum number of *m* values required. If a 
      value isn't passed, *mmax* is half the number of columns in ssphere
      minus one.

    Returns:
      ScalarCoefs: The object containing the coefficients of the scalar
      spherical harmonic transform.


    Raises:
      ValueError: If *nmax* and *mmax* are too large or *mmax* > *nmax*.

    """

    if nmax == None:
        nmax = ssphere.nrows - 2 
        mmax = int(ssphere.ncols / 2) - 1
    elif mmax == None:
        mmax = nmax

    if mmax > nmax:
        raise ValueError(err_msg['nmax_g_mmax'])

    if nmax >= ssphere.nrows - 1:
        raise ValueError(err_msg['nmax_too_lrg'])

    if mmax >= ssphere.ncols / 2:
        raise ValueError(err_msg['mmax_too_lrg'])

    dnrows = ssphere._dsphere.shape[0]
    ncols = ssphere._dsphere.shape[1]

    if np.mod(ncols, 2) == 1:
        raise ValueError(err_msg['ncols_even'])

    fdata = np.fft.fft2(ssphere._dsphere) / (dnrows * ncols)
    ops.fix_even_row_data_fc(fdata)
    
    fdata_extended = np.zeros([dnrows + 2, ncols], dtype=np.complex128)

    ops.pad_rows_fdata(fdata, fdata_extended)

    ops.sin_fc(fdata_extended)
    
    N = nmax + 1;
    NC = N + mmax * (2 * N - mmax - 1);
    sc = np.zeros(NC, dtype=np.complex128)
    # check if we are using c extended versions of the code or not
    if use_cext: 
        csphi.fc_to_sc(fdata_extended, sc, nmax, mmax)
    else:   
        sc = pysphi.fc_to_sc(fdata_extended, nmax, mmax)
                                
    return ScalarCoefs(sc, nmax, mmax)

def vspht(vsphere, nmax=None, mmax=None):
    """Returns a VectorCoefs object containt the vector spherical harmonic
    coefficients of the VectorPatternUniform object"""
    
    if nmax == None:
        nmax = vsphere.nrows - 2 
        mmax = int(vsphere.ncols / 2) - 1
    elif mmax == None:
        mmax = nmax

    if mmax > nmax:
        raise ValueError(err_msg['nmax_g_mmax'])

    if nmax >= vsphere.nrows - 1:
        raise ValueError(err_msg['nmax_too_lrg'])

    if mmax >= vsphere.ncols / 2:
        raise ValueError(err_msg['mmax_too_lrg'])

    dnrows = vsphere._tdsphere.shape[0]
    ncols = vsphere._tdsphere.shape[1]

    if np.mod(ncols, 2) == 1:
        raise ValueError(err_msg['ncols_even'])
        
    ft = np.fft.fft2(vsphere._tdsphere) / (dnrows * ncols)
    ops.fix_even_row_data_fc(ft)
    
    ft_extended = np.zeros([dnrows + 2, ncols], dtype=np.complex128)
    ops.pad_rows_fdata(ft, ft_extended)
    
    pt = np.fft.fft2(vsphere._pdsphere) / (dnrows * ncols)
    ops.fix_even_row_data_fc(pt)
    
    pt_extended = np.zeros([dnrows + 2, ncols], dtype=np.complex128)
    ops.pad_rows_fdata(pt, pt_extended)
    
    ftmp = np.copy(ft_extended)
    ptmp = np.copy(pt_extended)
    Lf1 = ops.sinLdot_fc(ft_extended, pt_extended)
    Lf2 = ops.sinLdot_fc(-1j * ptmp, 1j * ftmp)
    
    # check if we are using c extended versions of the code or not
    if use_cext: 
        N = nmax + 1;
        NC = N + mmax * (2 * N - mmax - 1);
        sc1 = np.zeros(NC, dtype=np.complex128)
        sc2 = np.zeros(NC, dtype=np.complex128)
        csphi.fc_to_sc(Lf1, sc1, nmax, mmax)
        csphi.fc_to_sc(Lf2, sc2, nmax, mmax)
    else:   
        sc1 = pysphi.fc_to_sc(Lf1, nmax, mmax)
        sc2 = pysphi.fc_to_sc(Lf2, nmax, mmax)

    vcoefs = VectorCoefs(sc1, sc2, nmax, mmax)

    nvec = np.zeros(nmax + 1, dtype=np.complex128)

    for n in xrange(1, nmax + 1):
        nvec[n] = 1.0 / np.sqrt(n * (n + 1.0))

    vcoefs.scoef1.window(nvec)
    vcoefs.scoef2.window(nvec)
        
    return vcoefs

def spht_slow(ssphere, nmax, mmax):
    """(PURE PYTHON) Transforms ScalarPatternUniform object *ssphere* 
    into a set of scalar spherical harmonics stored in ScalarCoefs.

    Example::

        >>> p = spherepy.random_patt_uniform(6, 8)
        >>> c = spherepy.spht(p)
        >>> spherepy.pretty_coefs(c)

    Args:
      ssphere (ScalarPatternUniform): The pattern to be transformed.

      nmax (int, optional): The maximum number of *n* values required. If a 
      value isn't passed, *nmax* is the number of rows in ssphere minus one.

      mmax (int, optional): The maximum number of *m* values required. If a 
      value isn't passed, *mmax* is half the number of columns in ssphere
      minus one.

    Returns:
      ScalarCoefs: The object containing the coefficients of the scalar
      spherical harmonic transform.


    Raises:
      ValueError: If *nmax* and *mmax* are too large or *mmax* > *nmax*.

    """

    if mmax > nmax:
        raise ValueError(err_msg['nmax_g_mmax'])

    nrows = ssphere._dsphere.shape[0]
    ncols = ssphere._dsphere.shape[1]

    if np.mod(nrows, 2) == 1 or np.mod(ncols, 2) == 1:
        raise ValueError(err_msg['ncols_even'])

    fdata = np.fft.fft2(ssphere._dsphere) / (nrows * ncols)
    ops.fix_even_row_data_fc(fdata)
    
    fdata_extended = np.zeros([nrows + 2, ncols], dtype=np.complex128)

    ops.pad_rows_fdata(fdata, fdata_extended)

    ops.sin_fc(fdata_extended)
    
    sc = pysphi.fc_to_sc(fdata_extended, nmax, mmax)
                                
    return ScalarCoefs(sc, nmax, mmax)

def ispht(scoefs, nrows=None, ncols=None):
    """Transforms ScalarCoefs object *scoefs* into a scalar pattern 
    ScalarPatternUniform.

    Example::

        >>> c = spherepy.random_coefs(3,3)
        >>> p = spherepy.ispht(c)
        >>> print(p)

    Args:
      scoefs (ScalarCoefs): The coefficients to be transformed to pattern
      space.

      nrows (int): The number of rows desired in the pattern.

      ncols (int): The number of columns desired in the pattern. This must be 
      an even number.

    Returns:
      ScalarPatternUniform: This is the pattern. It contains a NumPy array that
      can be viewed with *patt.cdata*.


    Raises:
      ValueError: Is raised if *ncols* isn't even.

      ValueError: Is raised if *nrows* < *nmax* + 2 or *ncols* < 2 * *mmax* + 2.

    """

    if nrows == None:
        nrows = scoefs.nmax + 2 

    if ncols == None:
        ncols = 2 * scoefs.mmax + 2

    if nrows <= scoefs.nmax:
        raise ValueError(err_msg['inverse_terr'])

    if ncols < 2 * scoefs.mmax + 2:
        raise ValueError(err_msg['inverse_terr'])

    dnrows = int(2 * nrows - 2)

    if np.mod(ncols, 2) == 1:
        raise ValueError(err_msg['ncols_even'])

    if use_cext: 
        fdata = np.zeros([dnrows, ncols], dtype=np.complex128)
        csphi.sc_to_fc(fdata, scoefs._vec, scoefs._nmax, scoefs._mmax)
    else:   
        fdata = pysphi.sc_to_fc(scoefs._vec,
                            scoefs._nmax,
                            scoefs._mmax,
                            dnrows, ncols)
    
    ds = np.fft.ifft2(fdata) * dnrows * ncols

    return ScalarPatternUniform(ds, doublesphere=True)

def vispht(vcoefs, nrows=None, ncols=None):



    if nrows == None:
        nrows = vcoefs.nmax + 2 

    if ncols == None:
        ncols = 2 * vcoefs.mmax + 2

    if nrows < vcoefs.nmax + 2:
        raise ValueError(err_msg['inverse_terr'])

    if ncols < 2 * vcoefs.mmax + 2:
        raise ValueError(err_msg['inverse_terr'])


    dnrows = int(2 * nrows - 2)

    if np.mod(ncols, 2) == 1:
        raise ValueError(err_msg['ncols_even'])

    nvec = np.zeros(vcoefs.nmax + 1, dtype=np.complex128)
    for n in xrange(1, vcoefs.nmax + 1):
        nvec[n] = 1.0 / np.sqrt(n * (n + 1.0))

    sc1 = vcoefs.scoef1.copy()
    sc2 = vcoefs.scoef2.copy()

    sc1.window(nvec)
    sc2.window(nvec)

    if use_cext: 
        fdata1 = np.zeros([dnrows, ncols], dtype=np.complex128)
        fdata2 = np.zeros([dnrows, ncols], dtype=np.complex128)
        csphi.sc_to_fc(fdata1, sc1._vec,
                       sc1.nmax, sc1.mmax)
        csphi.sc_to_fc(fdata2, sc2._vec,
                       sc1.nmax, sc1.mmax)
    else:   
        fdata1 = pysphi.sc_to_fc(sc1._vec,
                            sc1.nmax,
                            sc1.mmax,
                            dnrows, ncols)
        fdata2 = pysphi.sc_to_fc(sc2._vec,
                            sc2.nmax,
                            sc2.mmax,
                            dnrows, ncols)

    f1 = ops.L_fc(fdata1)
    f2 = ops.L_fc(fdata2)

    ftheta = f1[0] - 1j * f2[1]
    fphi = f1[1] + 1j * f2[0]
    
    dtheta = np.fft.ifft2(ftheta) * dnrows * ncols
    dphi = np.fft.ifft2(fphi) * dnrows * ncols

    return TransversePatternUniform(dtheta, dphi, doublesphere=True)

def ispht_slow(scoefs, nrows, ncols):
    """(PURE PYTHON) Transforms ScalarCoefs object *scoefs* into a scalar
    pattern ScalarPatternUniform.

    Example::

        >>> c = spherepy.random_coefs(3,3)
        >>> p = spherepy.ispht(c)
        >>> print(p)

    Args:
      scoefs (ScalarCoefs): The coefficients to be transformed to pattern
      space.

      nrows (int): The number of rows desired in the pattern.

      ncols (int): The number of columns desired in the pattern. This must be 
      an even number.

    Returns:
      ScalarPatternUniform: This is the pattern. It contains a NumPy array that
      can be viewed with *patt.cdata*.


    Raises:
      ValueError: Is raised if *ncols* isn't even.

    """

    dnrows = 2 * nrows - 2

    if np.mod(ncols, 2) == 1:
        raise ValueError(err_msg['ncols_even'])

    fdata = pysphi.sc_to_fc(scoefs._vec,
                            scoefs._nmax,
                            scoefs._mmax,
                            dnrows, ncols)
    
    ds = np.fft.ifft2(fdata) * dnrows * ncols

    return ScalarPatternUniform(ds, doublesphere=True)

def L2_coef(coef):

    if isinstance(coef, ScalarCoefs):
        return np.sqrt(np.sum(np.abs(coef._vec) ** 2))
    elif isinstance(coef, VectorCoefs):
        l1 = np.sqrt(np.sum(np.abs(coef.scoef1._vec) ** 2))
        l2 = np.sqrt(np.sum(np.abs(coef.scoef1._vec) ** 2))
        return np.sqrt(l1 ** 2 + l2 ** 2)
    
def L2_patt(patt):

    if isinstance(patt, ScalarPatternUniform):
        return np.sqrt(np.sum(np.abs(patt._dsphere) ** 2))
    elif isinstance(patt, TransversePatternUniform):
        l1 = np.sqrt(np.sum(np.abs(patt._tdsphere) ** 2))
        l2 = np.sqrt(np.sum(np.abs(patt._pdsphere) ** 2))
        return np.sqrt(l1 ** 2 + l2 ** 2)

def LInf_coef(coef):

    if isinstance(coef, ScalarCoefs):
        return np.max(np.abs(coef._vec))
    elif isinstance(coef, VectorCoefs):
        l1 = np.max(np.abs(coef.scoef1._vec))
        l2 = np.max(np.abs(coef.scoef2._vec))
        return np.sqrt(l1 ** 2 + l2 ** 2)
    
def LInf_patt(patt):

    if isinstance(patt, ScalarPatternUniform):
        return np.amax(np.abs(patt._dsphere))
    elif isinstance(patt, TransversePatternUniform):
        l1 = np.amax(np.abs(patt._tdsphere))
        l2 = np.amax(np.abs(patt._pdsphere))
        return np.sqrt(l1 ** 2 + l2 ** 2)
        

def sph_harmonic_tp(nrows, ncols, n, m):
    """Produces Ynm(theta,phi)

        theta runs from 0 to pi and has 'nrows' points.

        phi runs from 0 to 2*pi - 2*pi/ncols and has 'ncols' points.
    
    """
    nuvec = pysphi.ynunm(n, m, n + 1)
    num = nuvec[1:] * (-1) ** m
    num = num[::-1]
    yvec = np.concatenate([num, nuvec])
    
    theta = np.linspace(0.0, np.pi, nrows)
    phi = np.linspace(0.0, 2.0 * np.pi - 2.0 * np.pi / ncols)
    nu = np.array(range(-n, n + 1), dtype=np.complex128)
    
    out = np.zeros([nrows, ncols], dtype=np.complex128)

    for nn in xrange(0, nrows):
        exp_vec = np.exp(1j * nu * theta[nn])
        for mm in xrange(0, ncols):          
            vl = np.sum(yvec * exp_vec)
            out[nn, mm] = 1j ** m * np.exp(1j * m * phi[mm]) * vl

    return out





