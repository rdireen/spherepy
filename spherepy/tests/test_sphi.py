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


class TestSphi(TestCase):
    """Tests for the low level numeric functions"""
    
    def test_ynnm(self):
        """ynnm                
        """
        res = True
        for n in xrange(0,101):
            for m in xrange(-n,n):
                pyn = sp.pysphi.ynnm(n, m)
                cn = sp.csphi.ynnm(n, m)
                if (np.abs(pyn - cn)/np.abs(pyn) > 1e-13):
                    res = False    
        
        self.assertTrue(res)