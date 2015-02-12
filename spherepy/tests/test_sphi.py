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

from unittest import TestCase

import spherepy as sp
import numpy as np

Nmodes = 100

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
        for Nmax in xrange(10, 50):
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
        """spht"""
        
        Nrows = 303
        Ncols = 310
        
        
        
        
        self.assertTrue(True)
        
    def test_ispht(self):
        """ispht""" 
        
        Nmax = 100
        Nrows = 2000
        res = True
        
        for MM in xrange(1,5):
            Mmax = Nmax - MM
            c = sp.random_coefs(Nmax, Mmax)
            
            ppy = sp.ispht_slow(c, Nrows, Nrows)
            pc = sp.ispht_slow(c, Nrows, Nrows)
            
            rerr = np.sum(sp.abs(ppy - pc))/np.sum(sp.abs(ppy))
            
            if (rerr > 1e-13):
                res = False
        
        self.assertTrue(res)
        
   
        
        
                
                
        
        
        
        
        
        
        
