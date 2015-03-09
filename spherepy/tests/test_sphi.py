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
             test_sphi: Compare c extensions to python

Randy Direen
2/19/2015

These are tests that compare the c extensions to the python versions.

"""
from __future__ import division
from unittest import TestCase

import spherepy as sp
import numpy as np

# TODO: Change all xrange instances to range
# and do a 'from six.moves import range' here
from six.moves import xrange  # @UnresolvedImport

Nmodes = 20

class TestSphi(TestCase):
    """Tests for the low level numeric functions. This isn't a good test to 
    see if the numbers are right, but it is used for comparing the c code with
    the python code."""
    
    def test_ynnm(self):
        """::Compare ynnm within pysphi and csphi               
        """
        res = True
        for n in xrange(0, Nmodes + 1):
            for m in xrange(-n, n + 1):
                pyn = sp.pysphi.ynnm(n, m)
                cn = sp.csphi.ynnm(n, m)
                if (np.abs(pyn - cn) / np.abs(pyn) > 1e-13):
                    res = False    
        
        self.assertTrue(res)
        
    def test_ynunm(self):
        """::Compare ynunm within pysphi and csphi"""
        
        res = True
        for n in xrange(0, Nmodes + 1):
            for m in xrange(-n, n + 1):
                pyn = sp.pysphi.ynunm(n, m, n + 1)
                z = np.zeros(n + 1, dtype=np.float64)
                sp.csphi.ynunm(n, m, z)
                diff = np.sum(np.abs(pyn - z) ** 2)
                sm = np.sum(np.abs(pyn) ** 2)
                if (np.sqrt(diff / sm) > 1e-13):
                    res = False
                    
        self.assertTrue(res)
        
    def test_s_data(self):
        """::Compare s_data within pysphi and csphi"""
        
        res = True
        for Nmax in xrange(10, 20):
            for Nrows in xrange(2 * Nmax + 2, 20 * Nmax + 2, 2):
                Q = Nrows + Nmax + 1
                pyn = sp.pysphi.s_data(Nrows, Nmax, Q)
                z = np.zeros(Q, dtype=np.complex128)
                sp.csphi.SData(z, Nrows, Nmax)
                diff = np.sum(np.abs(pyn - z) ** 2)
                sm = np.sum(np.abs(pyn) ** 2)
                if (np.sqrt(diff / sm) > 1e-13):
                    res = False
        
        self.assertTrue(res)
        
    def test_spht(self):
        """::Compare spht within pysphi and csphi"""
        
        Nrows = 100
        Ncols = 200
        res = True
        
        for MM in xrange(2, 8, 2):
            NNcols = Ncols - MM
            p = sp.random_patt_uniform(Nrows, NNcols)
            
            spy = sp.spht_slow(p, 50, 46)
            sc = sp.spht(p, 50, 46)
            
            rerr = sp.compare_relative(spy, sc)
            if (rerr > 1e-13):
                res = False
        
        self.assertTrue(res)
        
    def test_ispht(self):
        """::Compare ispht within pysphi and csphi""" 
        
        Nmax = 50
        Nrows = 102
        res = True
        
        for MM in xrange(1, 3):
            Mmax = Nmax - MM
            c = sp.random_coefs(Nmax, Mmax)
            
            ppy = sp.ispht_slow(c, Nrows / 2, Nrows)
            pc = sp.ispht(c, Nrows / 2, Nrows)
            
            rerr = np.sum(sp.abs(ppy - pc)) / np.sum(sp.abs(ppy))
            
            if (rerr > 1e-13):
                res = False
                
        self.assertTrue(res)
        
   
        
        
                
                
        
        
        
        
        
        
        
