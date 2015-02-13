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
@newField description: Description
@newField revisions: Revisions
@newField departments: Departments
@newField applications: Applications

@author:
    rdireen

@organization:
    DireenTech Inc.

@description:

"""

#---------------------------------------------------------------------Built-ins
from functools import wraps
import numbers
import __init__

#---------------------------------------------------------------------3rd Party
import numpy as np

#------------------------------------------------------------------------Custom
import pysphi  # python versions of the low level routines
import csphi  # c extensions of the low level routines
import ops

#==============================================================================
# Global Declarations
#==============================================================================

scalar = 0
vector = 1

#==============================================================================
# Classes
#==============================================================================

class SpherePyError(Exception):
    pass

class ScalarCoefs(object):
    """Holds the scalar coefficients that represent a spherical pattern. The 
    function spht returns this object"""
    def __init__(self, vec, nmax, mmax):

        self._vec = vec
        self._nmax = nmax
        self._mmax = mmax

    @property
    def nmax(self):
        return self._nmax

    @property
    def mmax(self):
        return self._mmax

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
        return "scalar_coef((nmax = {0}, mmax = {1}) , {2} )".format(self.nmax,
                                                                     self.mmax,
                                                                     str(rs))

    def __str__(self):
        rs = self._reshape_m_vecs()
        st = ["scalar_coef(nmax = {0}, mmax = {1}):".format(self.nmax,
                                                          self.mmax)]
        for n, x in enumerate(rs):
            st.append("n = {0} :: {1}".format(n, x))
                      
        return '\n'.join(st)

    def __setitem__(self, arg, val):
        if len(arg) != 2:
            raise AttributeError("tuple must have 2 elements")

        if isinstance(arg[0], slice) and isinstance(arg[1], int):

            m = arg[1]

            if ((arg[0].start is None) and 
                (arg[0].step is None) and 
                (arg[0].stop is None)):

                if np.abs(m) > self._mmax:
                    raise AttributeError("m value out of bounds")

                idx_start = pysphi.mindx(m, self.nmax, self.mmax)
                idx_stop = idx_start + self.nmax - np.abs(m) + 1

                L = len(self._vec[idx_start:idx_stop])
                if len(val) != L:
                    raise AttributeError("dimension mismatch")
                else:
                    self._vec[idx_start:idx_stop] = val
            else:
                raise AttributeError("slicing operation not permitted")

        elif isinstance(arg[0], int) and isinstance(arg[1], slice):
            
            n = arg[0]

            if ((arg[1].start is None) and 
                (arg[1].step is None) and 
                (arg[1].stop is None)):

                if (n > self.nmax) or (n < 0):
                    raise AttributeError("n value out of bounds")

                if (len(val) != 2 * n + 1):
                    raise ArithmeticError("dimension mismatch")

                if n <= self.mmax:
                    nn = n
                else:
                    nn = self.mmax            
                for k, m in enumerate(xrange(-nn, nn + 1)):
                    self.__setitem__((n, m), val[k])

            else:
                raise AttributeError("slicing operation not permitted")

        elif isinstance(arg[0], int) and isinstance(arg[1], int):

            n = arg[0]
            m = arg[1]

            if (n < 0) or (n > self._nmax):
                raise AttributeError("n value out of bounds")

            if (np.abs(m) > n) or (np.abs(m) > self._mmax):
                raise AttributeError("m value out of bounds")

            idx = pysphi.mindx(m, self.nmax, self.mmax) + n - np.abs(m)

            self._vec[idx] = val

        else:
            raise AttributeError("what?, indexing method not recognized")

    def __getitem__(self, arg):

        if len(arg) != 2:
            raise AttributeError("tuple must have 2 elements")  

        if isinstance(arg[0], slice) and isinstance(arg[1], int):

            m = arg[1]

            if ((arg[0].start is None) and
               (arg[0].step is None) and
               (arg[0].stop is None)):

                if np.abs(m) > self._mmax:
                    raise AttributeError("m value out of bounds")

                idx_start = pysphi.mindx(m, self.nmax, self.mmax)
                idx_stop = idx_start + self.nmax - np.abs(m) + 1

                return self._vec[idx_start:idx_stop]
            else:
                raise AttributeError("slicing operation not permitted")

        elif isinstance(arg[0], int) and isinstance(arg[1], slice):
            
            n = arg[0]

            if ((arg[1].start is None) and
               (arg[1].step is None) and
               (arg[1].stop is None)):

                if (n > self.nmax) or (n < 0):
                    raise AttributeError("n value out of bounds")

                mlst = []
                if n <= self.mmax:
                    nn = n
                else:
                    nn = self.mmax            
                for m in xrange(-nn, nn + 1):
                    mlst.append(self.__getitem__((n, m)))

                return np.array(mlst, dtype=np.complex128)
            else:
                raise AttributeError("slicing operation not permitted")

        elif isinstance(arg[0], int) and isinstance(arg[1], int):

            n = arg[0]
            m = arg[1]

            if (n < 0) or (n > self._nmax):
                raise AttributeError("n value out of bounds")

            if (np.abs(m) > n) or (np.abs(m) > self._mmax):
                raise AttributeError("m value out of bounds")

            idx = pysphi.mindx(m, self.nmax, self.mmax) + n - np.abs(m)

            return self._vec[idx]

        else:
            raise AttributeError("what?, indexing method not recognized")

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
                    raise SpherePyError("ScalarCoefs could not be combined" + 
                                    " with sizes (%d,%d) and (%d,%d)" % \
                                    (self.nmax, self.mmax,
                                     scoef.nmax, scoef.mmax))
        
            elif isinstance(scoef, numbers.Number):
                return ScalarCoefs(func(self, self._vec, scoef), self.nmax,
                                   self.mmax)
            else:
                raise SpherePyError("cannot combine type with ScalarCoefs")
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
                raise SpherePyError("cannot add type to ScalarCoefs")
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
        return a / b

    @_scalar_coef_op_right
    def __rdiv__(self, a, b):
        return b / a  


class VectorCoefs(object):
    """Holds the vector coefficients that represent a vector spherical pattern. 
    The function vspht returns this object"""
    def __init__(self, vec1, vec2, nmax, mmax):
        
        if(vec1.shape != vec2.shape):
            raise SpherePyError("Shape of vec1 and vec2 must be the " + \
                                 "the same.")

        self.scoef1 = ScalarCoefs(vec1, nmax, mmax)
        self.scoef2 = ScalarCoefs(vec1, nmax, mmax)
        self._nmax = nmax
        self._mmax = mmax

    @property
    def nmax(self):
        return self._nmax

    @property
    def mmax(self):
        return self._mmax

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
        return "vector_coef(nmax = {0}, mmax = {1})".format(self.nmax,
                                                                     self.mmax)

    def __setitem__(self, arg, val):
        if len(arg) != 2:
            raise AttributeError("tuple must have 2 elements")
        
        if len(val) != 2:
            raise AttributeError("The value needs to be a tuple with " + \
                                 "length = 2 and containing two numbers.")

        if isinstance(arg[0], slice) and isinstance(arg[1], int):
            raise SpherePyError("Set scalar coefficients individually " + \
                                "(e.g. vcoef.scoef2[:,3] = vec)")

        elif isinstance(arg[0], int) and isinstance(arg[1], slice):
            raise SpherePyError("Set scalar coefficients individually " + \
                                "(e.g. vcoef.scoef2[5,:] = vec)")

        elif isinstance(arg[0], int) and isinstance(arg[1], int):
            n = arg[0]
            m = arg[1]

            try:
                self.scoef1[n, m] = val[0]
                self.scoef2[n, m] = val[1]
            except IndexError:
                SpherePyError("Set coefficients like this " + \
                              "vsc[n,m] = (val1,val2)")

        else:
            raise AttributeError("what?, indexing method not recognized")

    def __getitem__(self, arg):

        if len(arg) != 2:
            raise AttributeError("tuple must have 2 elements")  

        if isinstance(arg[0], slice) and isinstance(arg[1], int):

            m = arg[1]

            if ((arg[0].start is None) and
               (arg[0].step is None) and
               (arg[0].stop is None)):

                return (self.scoef1[:, m], self.scoef2[:, m])
            else:
                raise AttributeError("slicing operation not permitted")

        elif isinstance(arg[0], int) and isinstance(arg[1], slice):
            
            n = arg[0]

            if ((arg[1].start is None) and
               (arg[1].step is None) and
               (arg[1].stop is None)):
                
                return (self.scoef1[n, :], self.scoef2[n, :])
            
            else:
                raise AttributeError("slicing operation not permitted")

        elif isinstance(arg[0], int) and isinstance(arg[1], int):
            n = arg[0]
            m = arg[1]
            return (self.scoef1[n, m], self.scoef2[n, m])

        else:
            raise AttributeError("what?, indexing method not recognized")

    def _vector_coef_op_left(func):
        """decorator for operator overloading when VectorCoef is on the
        left"""
        @wraps(func)
        def verif(self, vcoef):
            if isinstance(vcoef, VectorCoefs):
                if len(self._vec) == len(vcoef._vec):
                    return VectorCoefs(func(self, self.scoef1._vec,
                                                  vcoef.scoef1._vec),
                                       func(self, self.scoef2._vec,
                                                  vcoef.scoef2._vec),
                                        self.nmax,
                                        self.mmax)
                else:
                    raise SpherePyError("VectorCoefs could not be combined" + 
                                    " with sizes (%d,%d) and (%d,%d)" % \
                                    (self.nmax, self.mmax,
                                     vcoef.nmax, vcoef.mmax))
        
            elif isinstance(vcoef, numbers.Number):
                return VectorCoefs(func(self, self.scoef1._vec, vcoef),
                                   func(self, self.scoef2._vec, vcoef),
                                   self.nmax,
                                   self.mmax)
            else:
                raise SpherePyError("cannot combine type with VectorCoefs")
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
                raise SpherePyError("cannot add type to VectorCoefs")
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
        return a / b

    @_vector_coef_op_right
    def __rdiv__(self, a, b):
        return b / a  


class ScalarPatternUniform(object):

    def __init__(self, cdata, doublesphere=False):

        if(doublesphere == False):
            self._dsphere = continue_sphere(cdata, 1)
        else:
            self._dsphere = cdata

    def __repr__(self):

        return self.array.__repr__()

    @property
    def nrows(self):
        return self._dsphere.shape[0] / 2 + 1 

    @property
    def ncols(self):
        return self._dsphere.shape[1]

    @property
    def shape(self):
        return (self.nrows, self.ncols)

    @property
    def doublesphere(self):
        return self._dsphere

    @property
    def array(self):
        return self._dsphere[0:self.nrows, :]

    @property
    def is_symmetric(self):
        """return true if the data is symmetric"""
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
        """decorator for operator overloading when ScalarPatternUniform is on 
        the left"""
        @wraps(func)
        def verif(self, patt):
            if isinstance(patt, ScalarPatternUniform):
                if self._dsphere.shape == patt._dsphere.shape:
                    return ScalarPatternUniform(func(self, self._dsphere,
                                                     patt._dsphere),
                                                doublesphere=True)
                else:
                    raise SpherePyError("ScalarPatternUniform could not be" + 
                                        " combined with sizes (%d,%d) and "
                                        + "(%d,%d)" % \
                                            (self.nrows, self.ncols,
                                            patt.nrows, patt.ncols))
        
            elif isinstance(patt, numbers.Number):
                return ScalarPatternUniform(func(self, self._dsphere, patt),
                                            doublesphere=True)
            else:
                raise SpherePyError("cannot combine type with" + 
                                   " ScalarPatternUniform")
        return verif

    def _scalar_pattern_uniform_op_right(func):
        """decorator for operator overloading when ScalarPatternUniform is on
        the right"""
        @wraps(func)
        def verif(self, patt):
            if isinstance(patt, numbers.Number):
                return ScalarPatternUniform(func(self, self._dsphere, patt),
                                   doublesphere=True)
            else:
                raise SpherePyError("cannot add type to ScalarPatternUniform")
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
        return a / b

    @_scalar_pattern_uniform_op_right
    def __rdiv__(self, a, b):
        return b / a  


class ScalarPatternNonUniform:
    pass

class VectorPatternUniform:
    
    def __init__(self, tcdata, pcdata, doublesphere=False):
        
        if(tcdata.shape != pcdata.shape):
            raise SpherePyError("Shape of tcdata and pcdata must be the " + \
                                 "the same")

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
        return self._tdsphere.shape[0] / 2 + 1 

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
                self._tdsphere[0:self.nrows, :])

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
            if isinstance(patt, VectorPatternUniform):
                if self._tdsphere.shape == patt._tdsphere.shape:
                    return VectorPatternUniform(func(self, self._tdsphere,
                                                     patt._tdsphere),
                                                func(self, self._pdsphere,
                                                     patt._pdsphere),
                                                doublesphere=True)
                else:
                    raise SpherePyError("VectorPatternUniform could not be" + 
                                        " combined with sizes (%d,%d) and "
                                        + "(%d,%d)" % \
                                            (self.nrows, self.ncols,
                                            patt.nrows, patt.ncols))
        
            elif isinstance(patt, numbers.Number):
                return VectorPatternUniform(func(self, self._tdsphere, patt),
                                            func(self, self._pdsphere, patt),
                                            doublesphere=True)
            else:
                raise SpherePyError("cannot combine type with" + 
                                   " VectorPatternUniform")
        return verif

    def _vector_pattern_uniform_op_right(func):
        """decorator for operator overloading when VectorPatternUniform is on
        the right"""
        @wraps(func)
        def verif(self, patt):
            if isinstance(patt, numbers.Number):
                return ScalarPatternUniform(func(self, self._tdsphere, patt),
                                            func(self, self._pdsphere, patt),
                                   doublesphere=True)
            else:
                raise SpherePyError("cannot add type to VectorPatternUniform")
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
        return a / b

    @_vector_pattern_uniform_op_right
    def __rdiv__(self, a, b):
        return b / a  


class VectorPatternNonUniform:
    pass

#==============================================================================
# Functions
#==============================================================================
        
#---------------------------------------------------------------------Factories
def zeros_coefs(nmax, mmax, coef_type=scalar):
    """return ScalarCoefs with all coefficients set to zeros"""

    if(mmax > nmax):
        raise SpherePyError("mmax must be less than nmax")

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
        raise SpherePyError("unknown coefficient type")

def ones_coefs(nmax, mmax, coef_type=scalar):
    """return ScalarCoefs with all coefficients set to ones"""

    if(mmax > nmax):
        raise SpherePyError("mmax must be less than nmax")
    
    if(coef_type == scalar):
        L = (nmax + 1) + mmax * (2 * nmax - mmax + 1)
        vec = np.ones(L, dtype=np.complex128)
        return  ScalarCoefs(vec, nmax, mmax)
    elif(coef_type == vector):
        L = (nmax + 1) + mmax * (2 * nmax - mmax + 1)
        vec1 = np.ones(L, dtype=np.complex128)
        vec2 = np.ones(L, dtype=np.complex128)
        return  VectorCoefs(vec1, vec2, nmax, mmax)
    else:
        raise SpherePyError("unknown coefficient type")

def random_coefs(nmax, mmax, mu=0.0, sigma=1.0, coef_type=scalar):
    """return ScalarCoefs with all coefficients random normal variables"""

    if(mmax > nmax):
        raise SpherePyError("mmax must be less than nmax")

    if(coef_type == scalar):
        L = (nmax + 1) + mmax * (2 * nmax - mmax + 1)
        vec = np.random.normal(mu, sigma, L) + \
              1j * np.random.normal(mu, sigma, L)
        return  ScalarCoefs(vec, nmax, mmax)
    elif(coef_type == vector):
        L = (nmax + 1) + mmax * (2 * nmax - mmax + 1)
        vec1 = np.random.normal(mu, sigma, L) + \
              1j * np.random.normal(mu, sigma, L)
        vec2 = np.random.normal(mu, sigma, L) + \
              1j * np.random.normal(mu, sigma, L)
        return  VectorCoefs(vec1, vec2, nmax, mmax)
    else:
        raise SpherePyError("unknown coefficient type")

def zeros_patt_uniform(nrows, ncols, patt_type=scalar):
    
    if np.mod(ncols, 2) == 1:
        raise SpherePyError("ncols must be even")

    if(patt_type == scalar):
        cdata = np.zeros((2 * nrows + 2, ncols), dtype=np.complex128)
        return ScalarPatternUniform(cdata, doublesphere=True)

    elif(patt_type == vector):
        tcdata = np.zeros((2 * nrows + 2, ncols), dtype=np.complex128)
        pcdata = np.zeros((2 * nrows + 2, ncols), dtype=np.complex128)
        return VectorPatternUniform(tcdata, pcdata, doublesphere=True)

    else:
        raise SpherePyError("unrecognized pattern type")

def ones_patt_uniform(nrows, ncols, patt_type=scalar):
    if np.mod(ncols, 2) == 1:
        raise SpherePyError("ncols must be even")

    if(patt_type == scalar):
        cdata = np.ones((nrows, ncols), dtype=np.complex128)
        return ScalarPatternUniform(cdata, doublesphere=False)

    elif(patt_type == vector):
        cdata = np.ones((nrows, ncols), dtype=np.complex128)
        return VectorPatternUniform(cdata, doublesphere=False)

    else:
        raise SpherePyError("unrecognized pattern type")

def random_patt_uniform(nrows, ncols, patt_type=scalar):

    if np.mod(ncols, 2) == 1:
        raise SpherePyError("ncols must be even")

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
        return VectorPatternUniform(vec1.reshape((nrows, ncols)),
                                    vec2.reshape((nrows, ncols)),
                                    doublesphere=False)

    else:
        raise SpherePyError("unrecognized pattern type")

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
    elif isinstance(patt, VectorPatternUniform):
        return patt.array
    else:
        raise SpherePyError("unrecognized type")
    
def abs(sobj):
    if isinstance(sobj, ScalarPatternUniform):
        return np.abs(sobj.array)
    elif isinstance(sobj, ScalarCoefs):
        return ScalarCoefs(np.abs(sobj._vec), sobj.nmax, sobj.mmax)
    elif isinstance(sobj, VectorPatternUniform):
        raise NotImplementedError()
    else:
        raise SpherePyError("unrecognized type")
    
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

def spht(ssphere, nmax, mmax):
    """Returns a ScalarCoefs object containing the spherical harmonic 
    coefficients of the ScalarPatternUniform object"""

    if mmax > nmax:
        raise SpherePyError("Mmax must be less than or equal to Nmax")

    nrows = ssphere._dsphere.shape[0]
    ncols = ssphere._dsphere.shape[1]

    if np.mod(nrows, 2) == 1 or np.mod(ncols, 2) == 1:
        raise SpherePyError("number of rows and columns in continued sphere" + 
                        " object must be even")

    fdata = np.fft.fft2(ssphere._dsphere) / (nrows * ncols)
    ops.fix_even_row_data_fc(fdata)
    
    fdata_extended = np.zeros([nrows + 2, ncols], dtype=np.complex128)

    ops.pad_rows_fdata(fdata, fdata_extended)

    ops.sin_fc(fdata_extended)
    
    N = nmax + 1;
    NC = N + mmax * (2 * N - mmax - 1);
    sc = np.zeros(NC, dtype=np.complex128)
    # check if we are using c extended versions of the code or not
    if __init__.use_cext: 
        csphi.fc_to_sc(fdata_extended, sc, nmax, mmax)
    else:   
        sc = pysphi.fc_to_sc(fdata_extended, nmax, mmax)
                                
    return ScalarCoefs(sc, nmax, mmax)

def vspht(vsphere, nmax, mmax):
    """Returns a VectorCoefs object containt the vector spherical harmonic
    coefficients of the VectorPatternUniform object"""
    
    if mmax > nmax:
        raise SpherePyError("Mmax must be less than or equal to Nmax")

    nrows = vsphere.scoef1_dsphere.shape[0]
    ncols = vsphere.scoef1_dsphere.shape[1]

    if np.mod(nrows, 2) == 1 or np.mod(ncols, 2) == 1:
        raise SpherePyError("number of rows and columns in continued sphere" + 
                        " object must be even")
        
    ft = np.fft.fft2(vsphere._tdsphere)
    ops.fix_even_row_data_fc(ft)
    
    ft_extended = np.zeros([nrows + 2, ncols], dtype=np.complex128)
    ops.pad_rows_fdata(ft, ft_extended)
    
    pt = np.fft.fft2(vsphere._pdsphere)
    ops.fix_even_row_data_fc(pt)
    
    pt_extended = np.zeros([nrows + 2, ncols], dtype=np.complex128)
    ops.pad_rows_fdata(pt, pt_extended)
    
    Lf1 = ops.sinLdot_fc(ft_extended, pt_extended)
    Lf2 = ops.sinLdot_fc(ft_extended, pt_extended)
    
    N = nmax + 1;
    NC = N + mmax * (2 * N - mmax - 1);
    sc1 = np.zeros(NC, dtype=np.complex128)
    sc2 = np.zeros(NC, dtype=np.complex128)
    # check if we are using c extended versions of the code or not
    if __init__.use_cext: 
        csphi.fc_to_sc(Lf1, sc1, nmax, mmax)
        csphi.fc_to_sc(Lf2, sc2, nmax, mmax)
    else:   
        sc1 = pysphi.fc_to_sc(Lf1, nmax, mmax)
        sc2 = pysphi.fc_to_sc(Lf1, nmax, mmax)
        
    return VectorCoefs(sc1, sc2, nmax, mmax)
    
    

def spht_slow(ssphere, nmax, mmax):
    """Returns a ScalarCoefs object containing the spherical harmonic 
    coefficients of the ssphere object"""

    if mmax > nmax:
        raise SpherePyError("Mmax must be less than or equal to Nmax")

    nrows = ssphere._dsphere.shape[0]
    ncols = ssphere._dsphere.shape[1]

    if np.mod(nrows, 2) == 1 or np.mod(ncols, 2) == 1:
        raise SpherePyError("number of rows and columns in continued sphere" + 
                        " object must be even")

    fdata = np.fft.fft2(ssphere._dsphere) / (nrows * ncols)
    ops.fix_even_row_data_fc(fdata)
    
    fdata_extended = np.zeros([nrows + 2, ncols], dtype=np.complex128)

    ops.pad_rows_fdata(fdata, fdata_extended)

    ops.sin_fc(fdata_extended)
    
    # check if we are using c extended versions of the code or not
    sc = pysphi.fc_to_sc(fdata_extended, nmax, mmax)
                                
    return ScalarCoefs(sc, nmax, mmax)

def ispht(scoefs, nrows, ncols):

    if np.mod(ncols, 2) == 1:
        raise SpherePyError("For consistency purposes, make sure Ncols " + 
                            "is even.")

    if __init__.use_cext: 
        fdata = np.zeros([nrows, ncols], dtype=np.complex128)
        csphi.sc_to_fc(fdata, scoefs._vec, scoefs._nmax, scoefs._mmax)
    else:   
        fdata = pysphi.sc_to_fc(scoefs._vec,
                            scoefs._nmax,
                            scoefs._mmax,
                            nrows, ncols)
    
    ds = np.fft.ifft2(fdata) * nrows * ncols

    return ScalarPatternUniform(ds, doublesphere=True)


def ispht_slow(scoefs, nrows, ncols):

    if np.mod(ncols, 2) == 1:
        raise SpherePyError("For consistency purposes, make sure Ncols " + 
                            "is even.")

    fdata = pysphi.sc_to_fc(scoefs._vec,
                            scoefs._nmax,
                            scoefs._mmax,
                            nrows, ncols)
    
    ds = np.fft.ifft2(fdata) * nrows * ncols

    return ScalarPatternUniform(ds, doublesphere=True)

def L2_coef(scoef):

    return np.sqrt(np.sum(np.abs(scoef._vec) ** 2))

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





