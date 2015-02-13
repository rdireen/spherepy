from unittest import TestCase

import spherepy as sp
import numpy as np
from spherepy.spherepy import SpherePyError


class TestScalarCoefsStructure(TestCase):
    
    def test_make_sure_mmax_not_greater_nmax(self):
        with self.assertRaises(ValueError):
            zz = sp.zeros_coefs(11, 13)
            
    
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
    
        print ds1.single_val
    
        print sp.array(ds1)
    
        zc = sp.zeros_coefs(10, 10)
        
        self.assertTrue(res)
        
class ScalarPatternUniform(TestCase):
        
    def test_make_sure_even_cols(self):
        """::Make sure I raise an error if I try to set ncol odd"""
        
        with self.assertRaises(SpherePyError):
            ww = sp.ones_coefs(12,13)
            