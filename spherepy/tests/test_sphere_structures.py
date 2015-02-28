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
Randy Direen
2/14/2015

I test the objects ScalarCoef, VectorCoef, ScalarPatternUniver, and 
VectorPatternUniform in this file. I test the exceptions, the indexing, and
how the overloaded operators work.

TODO: Test the single_val, doublesphere, continue, etc. functions in here as
well.

"""
from __future__ import division
from unittest import TestCase

import spherepy as sp
import numpy as np

#TODO: Change all xrange instances to range
#and do a 'from six.moves import range' here
from six.moves import xrange


class TestScalarCoefsStructure(TestCase):

    def test_pass_vec_wrong_size(self):
        """::raise an error if I try to pass a vec of wrong size"""
        vec = np.zeros(11, dtype=np.complex128)
        with self.assertRaises(ValueError):
            sc = sp.ScalarCoefs(vec,2,2)

    def test_create_sc(self):
        """::raise an error if I try to make mmax larger than nmax in
        ScalarCoefs"""
        vec = np.zeros(11, dtype=np.complex128)
        with self.assertRaises(ValueError):
            sc = sp.ScalarCoefs(vec,2,3)
    
    def test_make_sure_mmax_not_greater_nmax_z(self):
        """::raise an error if I try to set mmax > nmax: zeros"""  
        with self.assertRaises(ValueError):
            zz = sp.zeros_coefs(11, 13)
            
    def test_make_sure_mmax_not_greater_nmax_o(self):
        """::raise an error if I try to set mmax > nmax: ones""" 
        with self.assertRaises(ValueError):
            zz = sp.ones_coefs(11, 13)
            
    def test_make_sure_mmax_not_greater_nmax_r(self):
        """::raise an error if I try to set mmax > nmax: random""" 
        with self.assertRaises(ValueError):
            zz = sp.random_coefs(11, 13)
      
    def test_bounds_checking_access_m_greater_n(self):
        """::raise an error if I go out of bounds m > n access"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            a = zz[2,3]
            
    def test_bounds_checking_access_m_greater_mmax(self):
        """::raise an error if I go out of bounds m > mmax access"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            a = zz[11,-11]
            
    def test_bounds_checking_access_n_greater_nmax(self):
        """::raise an error if I go out of bounds n > nmax access"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            a = zz[12,-9]
     
    def test_two_many_indecies(self):
        """::see if something like this z[1,2,3]"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            zz[2,3,1] = 8
                   
    def test_bounds_checking_set_m_greater_n(self):
        """::raise an error if I go out of bounds m > n set"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            zz[2,3] = 8

    def test_bounds_checking_set_m_greater_mmax(self):
        """::raise an error if I go out of bounds m > mmax set"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            zz[11,11] = 8
            
    def test_bounds_checking_set_m_greater_mmax(self):
        """::raise an error if I go out of bounds m > mmax set"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            zz[11,11] = 8

    def test_bounds_checking_slice_n_not_right_size(self):
        """::raise an error I try to set slice to wrong size, zz[:,1] = vec TOO BIG"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            zz[:,1] = np.zeros(15, dtype = np.complex128)

    def test_bounds_checking_slice_n_slice_semi_m(self):
        """::raise an error, I can't slice like that yet z[:,-1:1] = vec"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            zz[:,1:2] = np.zeros(15, dtype = np.complex128)

    def test_indexing_not_recognized(self):
        """::raise an error, z[apple,:]"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            zz[np,:] = 8

    def test_bounds_checking_set_n_greater_nmax(self):
        """::raise an error if I go out of bounds n > nmax set"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            zz[12,10] = 8

    def test_bounds_checking_set_larger_than_one_item(self):
        """::raise an error if I set with larger tuple"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            zz[12,10] = (8, 0)

    def test_slice_to_sc1(self):
        """::raise an error if I slice incorrectly1"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            sc = zz[0:5,10]

    def test_too_many_idx(self):
        """::raise an error too many indexes"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            sc = zz[0:5,10,:]

    def test_slice_to_sc2(self):
        """::raise an error if I slice incorrectly2"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            sc = zz[1:5,:]
            
    def test_add_diff_sizes(self):
        """::can't add two ScalarCoefs that are different """
        z1 = sp.random_coefs(11, 10)
        z2 = sp.random_coefs(12, 10)
        with self.assertRaises(ValueError):
            a = z1 + z2
            
    def test_sub_diff_sizes(self):
        """::can't sub two ScalarCoefs that are different """
        z1 = sp.random_coefs(11, 10)
        z2 = sp.random_coefs(12, 10)
        with self.assertRaises(ValueError):
            a = z1 - z2
            
    def test_mult_diff_sizes(self):
        """::can't mult two ScalarCoefs that are different """
        z1 = sp.random_coefs(11, 10)
        z2 = sp.random_coefs(12, 10)
        with self.assertRaises(ValueError):
            a = z1 * z2
            
    def test_div_diff_sizes(self):
        """::can't div two ScalarCoefs that are different """
        z1 = sp.random_coefs(11, 10)
        z2 = sp.random_coefs(12, 10)
        with self.assertRaises(ValueError):
            a = z1 / z2

    def test_bounds_checking_n_slice_m(self):
        """::exercise pretty_coefs"""
        zz = sp.random_coefs(11, 10)
        a = sp.pretty_coefs(zz)
        self.assertTrue(True)

    def test_slice_ms(self):
        """::this should work a = z[1,:]"""
        zz = sp.random_coefs(11, 10)
        a = zz[1,:]
        self.assertTrue(True)

    def test_both_sliced(self):
        """::this should work a = z[0:3,:]"""
        zz = sp.random_coefs(11, 10)
        a = zz[0:3,:]
        self.assertTrue(True)

    def test_size(self):
        """::make sure c.size works."""
        z1 = sp.random_coefs(11, 10)
        N = z1.nmax + 1;
        NC = N + z1.mmax * (2 * N - z1.mmax - 1);

        res = False
        if z1.size == NC:
            res = True

        self.assertTrue(res)

    def test___array_2d_repr(self):
        """::excercise _array_2d_repr"""
        z1 = sp.random_coefs(11, 10)
        x = z1._array_2d_repr()
        self.assertTrue(True)

    def test___reshape_n_vecs(self):
        """::excercise _reshape_n_vecs"""
        z1 = sp.random_coefs(11, 10)
        nv = z1._reshape_n_vecs()

        res = False
        if len(nv) == 2*z1.mmax + 1:
            res = True

        self.assertTrue(res)

    def test___reshape_m_vecs(self):
        """::excercise _reshape_m_vecs"""
        z1 = sp.random_coefs(11, 10)
        nv = z1._reshape_m_vecs()

        res = False
        if len(nv) == z1.nmax + 1:
            res = True

        self.assertTrue(res)

    def test__repr_str(self):
        """::excercise __repr__ and __str__"""

        z1 = sp.random_coefs(2, 1)
        a = z1.__repr__()
        a = str(z1)

        self.assertTrue(True)

        
                 
            
    def test_add_same_sizes(self):
        """::can add two ScalarCoefs that are same sized """
        z1 = sp.random_coefs(11, 10)
        z2 = sp.random_coefs(11, 10)
        a = z1 + z2
        b = z1 + z2
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
        
            
    def test_sub_same_sizes(self):
        """::can sub two ScalarCoefs that are same sized """
        z1 = sp.random_coefs(11, 10)
        z2 = sp.random_coefs(11, 10)
        a = z1 - z2
        b = z1 - z2
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_mult_same_sizes(self):
        """::can mult two ScalarCoefs that are same sized """
        z1 = sp.random_coefs(11, 10)
        z2 = sp.random_coefs(11, 10)
        a = z1 * z2
        b = z1 * z2
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_div_same_sizes(self):
        """::can add two ScalarCoefs that are same sized """
        z1 = sp.random_coefs(11, 10)
        z2 = sp.random_coefs(11, 10)
        a = z1 / (z2 + 0.0001)
        b = z1 / (z2 + 0.0001)
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
        
    def test_add_constant_left_right(self):
        """::can add a scalar to ScalarCoefs from both sides"""
        z1 = sp.random_coefs(11, 10)
        
        a = z1 + 1j*1.1
        b = z1 + 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        a = 1j*1.1 + z1
        b = 1j*1.1 + z1
        c = a - b
        
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
        
            
    def test_sub_constant_left_right(self):
        """::can sub a scalar to ScalarCoefs from both sides"""
        z1 = sp.random_coefs(11, 10)
        
        a = z1 - 1j*1.1
        b = z1 - 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        a = 1j*1.1 - z1
        b = 1j*1.1 - z1
        c = a - b
        
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_mult_constant_left_right(self):
        """::can mult a scalar to ScalarCoefs from both sides"""
        z1 = sp.random_coefs(11, 10)
        
        a = z1 * 1j*1.1
        b = z1 * 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        a = 1j*1.1 * z1
        b = 1j*1.1 * z1
        c = a - b
        
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_div_constant_left_right(self):
        """::can mult a scalar to ScalarCoefs from both sides"""
        z1 = sp.random_coefs(11, 10)
        
        a = z1 / 1j*1.1
        b = z1 / 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        a = 1j*1.1 / z1
        b = 1j*1.1 / z1
        c = a - b
        
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_scalar_coefs_ops(self):
        """::Check to make sure the scalar coefficient ops are working            
        """
        
        ww = sp.ones_coefs(12,8)
        rr = sp.random_coefs(10, 10) + 1j * sp.random_coefs(10, 10)
        zz = sp.zeros_coefs(21,18)

        ss = 1j - 2 / (sp.zeros_coefs(5, 4) + .001) + \
                4 * sp.random_coefs(5, 4) / 6.0
        ss += sp.ones_coefs(5, 4, coef_type=sp.scalar)
    
        qq = 1 + 2 * sp.ones_coefs(3, 3) / 4 * sp.ones_coefs(3, 3) - 2
    
        ds1 = 1j + 4 * sp.ones_patt_uniform(5, 8) / 3 - \
            sp.ones_patt_uniform(5, 8) - 1 + \
            10 * sp.random_patt_uniform(5, 8)
    
        ds1.single_val
    
        sp.array(ds1)
    
        zc = sp.zeros_coefs(10, 10);
        
        self.assertTrue(True)
        
class TestVectorCoefsStructure(TestCase):

    def test_create_vc(self):
        """::raise an error if I try to pass a vec1 and vec2 of wrong size"""
        vec = np.zeros(11, dtype=np.complex128)
        with self.assertRaises(ValueError):
            vc = sp.VectorCoefs(vec,vec,2,2)
    
    def test_make_sure_mmax_not_greater_nmax_z(self):
        """::raise an error if I try to set mmax > nmax: zeros"""  
        with self.assertRaises(ValueError):
            zz = sp.zeros_coefs(11, 13,coef_type=sp.vector)
            
    def test_make_sure_mmax_not_greater_nmax_o(self):
        """::raise an error if I try to set mmax > nmax: ones""" 
        with self.assertRaises(ValueError):
            zz = sp.ones_coefs(11, 13,coef_type=sp.vector)
            
    def test_make_sure_mmax_not_greater_nmax_r(self):
        """::raise an error if I try to set mmax > nmax: random""" 
        with self.assertRaises(ValueError):
            zz = sp.random_coefs(11, 13)
     
    def test_monopole_access(self):
        """::raise an error if I try to access monopole"""
        zz = sp.random_coefs(11, 10,coef_type=sp.vector)
        with self.assertRaises(AttributeError):
            a = zz[0,0]

    def test_monopole_set(self):
        """::raise an error if I try to set monopole"""
        zz = sp.random_coefs(11, 10,coef_type=sp.vector)
        with self.assertRaises(AttributeError):
            zz[0,0] = (1,2)
             
    def test_bounds_checking_access_m_greater_n(self):
        """::raise an error if I go out of bounds m > n access"""
        zz = sp.random_coefs(11, 10,coef_type=sp.vector)
        with self.assertRaises(AttributeError):
            a = zz[2,3]
            
    def test_bounds_checking_access_m_greater_mmax(self):
        """::raise an error if I go out of bounds m > mmax access"""
        zz = sp.random_coefs(11, 10,coef_type=sp.vector)
        with self.assertRaises(AttributeError):
            a = zz[11,-11]
            
    def test_bounds_checking_access_n_greater_nmax(self):
        """::raise an error if I go out of bounds n > nmax access"""
        zz = sp.random_coefs(11, 10,coef_type=sp.vector)
        with self.assertRaises(AttributeError):
            a = zz[12,-9]
            
    def test_bounds_checking_set_m_greater_n(self):
        """::raise an error if I go out of bounds m > n set"""
        zz = sp.random_coefs(11, 10,coef_type=sp.vector)
        with self.assertRaises(AttributeError):
            zz[2,3] = 8
            
    def test_bounds_checking_set_m_greater_mmax(self):
        """::raise an error if I go out of bounds m > mmax set"""
        zz = sp.random_coefs(11, 10,coef_type=sp.vector)
        with self.assertRaises(AttributeError):
            zz[11,11] = (8,1j)
            
    def test_bounds_checking_set_n_greater_nmax(self):
        """::raise an error if I go out of bounds n > nmax set"""
        zz = sp.random_coefs(11, 10,coef_type=sp.vector)
        with self.assertRaises(AttributeError):
            zz[12,10] = (8,1j)

    def test_add_diff_sizes(self):
        """::can't add two VectorCoefs that are different """
        z1 = sp.random_coefs(11, 10,coef_type=sp.vector)
        z2 = sp.random_coefs(12, 10,coef_type=sp.vector)
        with self.assertRaises(ValueError):
            a = z1 + z2
            
    def test_sub_diff_sizes(self):
        """::can't sub two VectorCoefs that are different """
        z1 = sp.random_coefs(11, 10,coef_type=sp.vector)
        z2 = sp.random_coefs(12, 10,coef_type=sp.vector)
        with self.assertRaises(ValueError):
            a = z1 - z2
            
    def test_mult_diff_sizes(self):
        """::can't mult two ScalarCoefs that are different """
        z1 = sp.random_coefs(11, 10,coef_type=sp.vector)
        z2 = sp.random_coefs(12, 10,coef_type=sp.vector)
        with self.assertRaises(ValueError):
            a = z1 * z2
            
    def test_div_diff_sizes(self):
        """::can't div two VectorCoefs that are different """
        z1 = sp.random_coefs(11, 10,coef_type=sp.vector)
        z2 = sp.random_coefs(12, 10,coef_type=sp.vector)
        with self.assertRaises(ValueError):
            a = z1 / z2

    def test_size(self):
        """::can add get the size? """
        z1 = sp.random_coefs(11, 10,coef_type=sp.vector)
        z1.size    
        self.assertTrue(True)
            
    def test_add_same_sizes(self):
        """::can add two VectorCoefs that are same sized """
        z1 = sp.random_coefs(11, 10,coef_type=sp.vector)
        z2 = sp.random_coefs(11, 10,coef_type=sp.vector)
        a = z1 + z2
        b = z1 + z2
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
        
            
    def test_sub_same_sizes(self):
        """::can sub two VectorCoefs that are same sized """
        z1 = sp.random_coefs(11, 10,coef_type=sp.vector)
        z2 = sp.random_coefs(11, 10,coef_type=sp.vector)
        a = z1 - z2
        b = z1 - z2
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_mult_same_sizes(self):
        """::can mult two VectorCoefs that are same sized """
        z1 = sp.random_coefs(11, 10,coef_type=sp.vector)
        z2 = sp.random_coefs(11, 10,coef_type=sp.vector)
        a = z1 * z2
        b = z1 * z2
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_div_same_sizes(self):
        """::can add two VectorCoefs that are same sized """
        z1 = sp.random_coefs(11, 10,coef_type=sp.vector)
        z2 = sp.random_coefs(11, 10,coef_type=sp.vector)
        a = z1 / (z2 + 0.0001)
        b = z1 / (z2 + 0.0001)
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
        
    def test_add_constant_left_right(self):
        """::can add a scalar to VectorCoefs from both sides"""
        z1 = sp.random_coefs(11, 10,coef_type=sp.vector)
        
        a = z1 + 1j*1.1
        b = z1 + 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        a = 1j*1.1 + z1
        b = 1j*1.1 + z1
        c = a - b
        
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
        
            
    def test_sub_constant_left_right(self):
        """::can sub a scalar to VectorCoefs from both sides"""
        z1 = sp.random_coefs(11, 10,coef_type=sp.vector)
        
        a = z1 - 1j*1.1
        b = z1 - 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        a = 1j*1.1 - z1
        b = 1j*1.1 - z1
        c = a - b
        
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_mult_constant_left_right(self):
        """::can mult a scalar to VectorCoefs from both sides"""
        z1 = sp.ones_coefs(11, 10,coef_type=sp.vector)
        
        a = z1 * 1j*1.1
        b = z1 * 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        a = 1j*1.1 * z1
        b = 1j*1.1 * z1
        c = a - b
        
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_div_constant_left_right(self):
        """::can mult a scalar to VectorCoefs from both sides"""
        z1 = sp.random_coefs(11, 10,coef_type=sp.vector)
        
        a = z1 / 1j*1.1
        b = z1 / 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        a = 1j*1.1 / z1
        b = 1j*1.1 / z1
        c = a - b
        
        if sp.L2_coef(c) / sp.L2_coef(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
   

class TestScalarPatternUniform(TestCase):
        
    def test_make_sure_even_cols_z(self):
        """::raise an error if I try to set ncol odd zeros"""       
        with self.assertRaises(ValueError):
            ww = sp.zeros_patt_uniform(13, 17)
            
    def test_make_sure_even_cols_o(self):
        """::raise an error if I try to set ncol odd ones"""       
        with self.assertRaises(ValueError):
            ww = sp.ones_patt_uniform(13, 17)
            
    def test_make_sure_even_cols_r(self):
        """::raise an error if I try to set ncol odd random"""       
        with self.assertRaises(ValueError):
            ww = sp.random_patt_uniform(13, 17)
            
    def test_add_diff_sizes(self):
        """::can't add two ScalarPatternUniform that are different """
        z1 = sp.zeros_patt_uniform(11, 10)
        z2 = sp.random_patt_uniform(12, 10)
        with self.assertRaises(ValueError):
            a = z1 + z2
            
    def test_sub_diff_sizes(self):
        """::can't sub two ScalarPatternUniform that are different """
        z1 = sp.zeros_patt_uniform(11, 10)
        z2 = sp.random_patt_uniform(11, 12)
        with self.assertRaises(ValueError):
            a = z1 - z2
            
    def test_mult_diff_sizes(self):
        """::can't mult two ScalarPatternUniform that are different """
        z1 = sp.ones_patt_uniform(12, 12)
        z2 = sp.random_patt_uniform(12, 10)
        with self.assertRaises(ValueError):
            a = z1 * z2
            
    def test_div_diff_sizes(self):
        """::can't div two ScalarPatternUniform that are different """
        z1 = sp.zeros_patt_uniform(11, 10)
        z2 = sp.ones_patt_uniform(12, 12)
        with self.assertRaises(ValueError):
            a = z1 / z2


    def test_single_val(self):
        """::test single_val method"""
        z2 = sp.random_patt_uniform(11, 10)
        if z2.single_val > 1e-13:
            self.assertTrue(False)

        self.assertTrue(True)

    def test_is_symmetric(self):
        """::test single_val method"""
        z2 = sp.random_patt_uniform(11, 10)
        self.assertTrue(z2.is_symmetric)


    def test_array_doublesphere(self):
        """::exercise array and doublesphere"""
        z2 = sp.random_patt_uniform(11, 10)
        a = z2.array
        d = z2.doublesphere
        s = z2.shape
        nr = z2.nrows
        nc = z2.ncols
        self.assertTrue(True)

    
            
    def test_add_same_sizes(self):
        """::can add two ScalarPatternUniform that are same sized """
        z1 = sp.zeros_patt_uniform(11, 10)
        z2 = sp.random_patt_uniform(11, 10)
        a = z1 + z2
        b = z1 + z2
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
        
            
    def test_sub_same_sizes(self):
        """::can sub two ScalarPatternUniform that are same sized """
        z1 = sp.zeros_patt_uniform(11, 10)
        z2 = sp.random_patt_uniform(11, 10)
        a = z1 - z2
        b = z1 - z2
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_mult_same_sizes(self):
        """::can mult two ScalarPatternUniform that are same sized """
        z1 = sp.ones_patt_uniform(11, 10)
        z2 = sp.random_patt_uniform(11, 10)
        a = z1 * z2
        b = z1 * z2
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_div_same_sizes(self):
        """::can add two ScalarPatternUniform that are same sized """
        z1 = sp.ones_patt_uniform(11, 10)
        z2 = sp.random_patt_uniform(11, 10)
        a = z1 / (z2 + 0.0001)
        b = z1 / (z2 + 0.0001)
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
        
    def test_add_constant_left_right(self):
        """::can add a scalar to ScalarPatternUniform from both sides"""
        z1 = sp.ones_patt_uniform(11, 10)
        
        a = z1 + 1j*1.1
        b = z1 + 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        a = 1j*1.1 + z1
        b = 1j*1.1 + z1
        c = a - b
        
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
        
            
    def test_sub_constant_left_right(self):
        """::can sub a scalar to ScalarPatternUniform from both sides"""
        z1 = sp.ones_patt_uniform(11, 10)
        
        a = z1 - 1j*1.1
        b = z1 - 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        a = 1j*1.1 - z1
        b = 1j*1.1 - z1
        c = a - b
        
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_mult_constant_left_right(self):
        """::can mult a scalar to ScalarPatternUniform from both sides"""
        z1 = sp.ones_patt_uniform(11, 10)
        
        a = z1 * 1j*1.1
        b = z1 * 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        a = 1j*1.1 * z1
        b = 1j*1.1 * z1
        c = a - b
        
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_div_constant_left_right(self):
        """::can mult a scalar to ScalarPatternUniform from both sides"""
        z1 = sp.ones_patt_uniform(11, 10)
        
        a = z1 / 1j*1.1
        b = z1 / 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        a = 1j*1.1 / z1
        b = 1j*1.1 / z1
        c = a - b
        
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
        
class TestVectorPatternUniform(TestCase):
        
    def test_make_sure_even_cols_z(self):
        """::raise an error if I try to set ncol odd zeros"""       
        with self.assertRaises(ValueError):
            ww = sp.zeros_patt_uniform(13, 17, patt_type = sp.vector)
            
    def test_make_sure_even_cols_o(self):
        """::raise an error if I try to set ncol odd ones"""       
        with self.assertRaises(ValueError):
            ww = sp.ones_patt_uniform(13, 17, patt_type = sp.vector)
            
    def test_make_sure_even_cols_r(self):
        """::raise an error if I try to set ncol odd random"""       
        with self.assertRaises(ValueError):
            ww = sp.random_patt_uniform(13, 17, patt_type = sp.vector)
            
    def test_add_diff_sizes(self):
        """::can't add two VectorPatternUniform that are different """
        z1 = sp.zeros_patt_uniform(11, 10, patt_type = sp.vector)
        z2 = sp.random_patt_uniform(12, 10, patt_type = sp.vector)
        with self.assertRaises(ValueError):
            a = z1 + z2
            
    def test_sub_diff_sizes(self):
        """::can't sub two VectorPatternUniform that are different """
        z1 = sp.zeros_patt_uniform(11, 10, patt_type = sp.vector)
        z2 = sp.random_patt_uniform(11, 12, patt_type = sp.vector)
        with self.assertRaises(ValueError):
            a = z1 - z2
            
    def test_mult_diff_sizes(self):
        """::can't mult two VectorPatternUniform that are different """
        z1 = sp.ones_patt_uniform(12, 12, patt_type = sp.vector)
        z2 = sp.random_patt_uniform(12, 10, patt_type = sp.vector)
        with self.assertRaises(ValueError):
            a = z1 * z2
            
    def test_div_diff_sizes(self):
        """::can't div two VectorPatternUniform that are different """
        z1 = sp.zeros_patt_uniform(11, 10, patt_type = sp.vector)
        z2 = sp.ones_patt_uniform(12, 12, patt_type = sp.vector)
        with self.assertRaises(ValueError):
            a = z1 / z2

    def test_array_doublesphere(self):
        """::exercise array and doublesphere"""
        z1 = sp.random_patt_uniform(11, 10, patt_type = sp.vector)
        a = z1.array
        d = z1.doublesphere
        s = z1.shape
        nr = z1.nrows
        nc = z1.ncols
        self.assertTrue(True)


    #TODO: This needs to be addressed
    def test_single_val(self):
        """::test single_val method"""
        z1 = sp.random_patt_uniform(11, 10, patt_type = sp.vector)
        a = z1.single_val
        #if z1.single_val > 1e-10:
         #   self.assertTrue(False)

        self.assertTrue(True)

    def test_is_symmetric(self):
        """::test single_val method"""
        z2 = sp.random_patt_uniform(11, 10, patt_type = sp.vector)
        self.assertTrue(z2.is_symmetric)
            
    def test_add_same_sizes(self):
        """::can add two VectorPatternUniform that are same sized """
        z1 = sp.zeros_patt_uniform(11, 10, patt_type = sp.vector)
        z2 = sp.random_patt_uniform(11, 10, patt_type = sp.vector)
        a = z1 + z2
        b = z1 + z2
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
        
            
    def test_sub_same_sizes(self):
        """::can sub two VectorPatternUniform that are same sized """
        z1 = sp.ones_patt_uniform(11, 10, patt_type = sp.vector)
        z2 = sp.random_patt_uniform(11, 10, patt_type = sp.vector)
        a = z1 - z2
        b = z1 - z2
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_mult_same_sizes(self):
        """::can mult two VectorPatternUniform that are same sized """
        z1 = sp.ones_patt_uniform(11, 10, patt_type = sp.vector)
        z2 = sp.random_patt_uniform(11, 10, patt_type = sp.vector)
        a = z1 * z2
        b = z1 * z2
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_div_same_sizes(self):
        """::can add two VectorPatternUniform that are same sized """
        z1 = sp.ones_patt_uniform(11, 10, patt_type = sp.vector)
        z2 = sp.random_patt_uniform(11, 10, patt_type = sp.vector)
        a = z1 / (z2 + 0.0001)
        b = z1 / (z2 + 0.0001)
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
        
    def test_add_constant_left_right(self):
        """::can add a scalar to VectorPatternUniform from both sides"""
        z1 = sp.ones_patt_uniform(11, 10, patt_type = sp.vector)
        
        a = z1 + 1j*1.1
        b = z1 + 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        a = 1j*1.1 + z1
        b = 1j*1.1 + z1
        c = a - b
        
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
        
            
    def test_sub_constant_left_right(self):
        """::can sub a scalar to VectorPatternUniform from both sides"""
        z1 = sp.ones_patt_uniform(11, 10, patt_type = sp.vector)
        
        a = z1 - 1j*1.1
        b = z1 - 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        a = 1j*1.1 - z1
        b = 1j*1.1 - z1
        c = a - b
        
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_mult_constant_left_right(self):
        """::can mult a scalar to VectorPatternUniform from both sides"""
        z1 = sp.random_patt_uniform(11, 10, patt_type = sp.vector)
        
        a = z1 * 1j*1.1
        b = z1 * 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        a = 1j*1.1 * z1
        b = 1j*1.1 * z1
        c = a - b
        
        
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
    def test_div_constant_left_right(self):
        """::can div a scalar to VectorPatternUniform from both sides"""
        z1 = sp.random_patt_uniform(11, 10, patt_type = sp.vector)
        
        a = z1 / 1j*1.1
        b = z1 / 1j*1.1
        c = a - b
        
        res = True
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        a = 1j*1.1 / z1
        b = 1j*1.1 / z1
        c = a - b
        
        if sp.L2_patt(c) / sp.L2_patt(b) > 1e-13:
            res = False
            
        self.assertTrue(res)
            
class TestMiscRoutines(TestCase):
     
    def test_misc(self):
        """::exercies misc functions"""

        z1 = sp.random_patt_uniform(11, 10, patt_type = sp.scalar)

        a = sp.abs(z1)
        a = sp.mag(z1)
        a = sp.mag2(z1)
        a = sp.L2_patt(z1)
        a = sp.LInf_patt(z1)

        z1 = sp.random_patt_uniform(11, 10, patt_type = sp.vector)
        
        a = sp.mag(z1)
        a = sp.mag2(z1)
        a = sp.L2_patt(z1)
        a = sp.LInf_patt(z1)

        z1 = sp.random_coefs(11, 10, coef_type = sp.scalar)
        
        a = sp.abs(z1)
        a = sp.mag(z1)
        a = sp.mag2(z1)
        a = sp.L2_coef(z1)
        a = sp.LInf_coef(z1)

        z1 = sp.random_coefs(11, 10, coef_type = sp.vector)
        
        a = sp.mag(z1)
        a = sp.mag2(z1)
        a = sp.L2_coef(z1)
        a = sp.LInf_coef(z1)
            
        self.assertTrue(True)   

    def test_misc_raise1(self):
        """::exercies raise TypError abs vpatt"""
        z1 = sp.random_patt_uniform(11, 10, patt_type = sp.vector)
        with self.assertRaises(TypeError):
            a = sp.abs(z1)

    def test_misc_raise1(self):
        """::exercies misc functions"""   
        z1 = sp.random_coefs(11, 10, coef_type = sp.vector)
        with self.assertRaises(TypeError):
            a = sp.abs(z1)
      
            
    
            