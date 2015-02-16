from unittest import TestCase

import spherepy as sp
import numpy as np

"""
Randy Direen
2/14/2015

I test the objects ScalarCoef, VectorCoef, ScalarPatternUniver, and 
VectorPatternUniform in this file. I test the exceptions, the indexing, and
how the overloaded operators work.

TODO: Test the single_val, doublesphere, continue, etc. functions in here as
well.

"""


class TestScalarCoefsStructure(TestCase):
    
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
            
    def test_bounds_checking_set_n_greater_nmax(self):
        """::raise an error if I go out of bounds n > nmax set"""
        zz = sp.random_coefs(11, 10)
        with self.assertRaises(AttributeError):
            zz[12,10] = 8
            
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
        z1 = sp.random_coefs(11, 10,coef_type=sp.vector)
        
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
            
    
    
            