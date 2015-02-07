from unittest import TestCase

import spherepy as sp
import numpy as np



class TestSBessel(TestCase):
    """Tests for the Spherical Bessel functions"""
    
    def test_sbesselj_sums(self):
        """::Test the jn(z) functions using the summation formula:
        
            Inf
            sum  (2*n+1) * jn(z)**2 = 1
            n=0
            
        I test arguments from around 1.0 to 10,000 to make sure i get the sum 
        to within machine precision.                 
        """
        
        res = True
        try:
            for z in np.linspace(1.0, 10000.0,10):
                rnd = np.random.normal(0,100)
                ar = np.abs(rnd) + 1.1
                s = sp.sbesselj_sum(ar,int(np.floor(z+400)))
                if s > 1e-13:
                    res = False
        except:
            res = False
        
        self.assertTrue(res)


    def test_cross_product(self):
        """::Uses the cross-product relationship to test the routines: 

            j[n+1]y[n] - j[n]y[n+1] = 1 / (z**2)
    
        where j and y are sbesselj or sbessely vectors for a particular z.
        Doing this provides a check to see if the bessel functions are being
        calculated correctly.
    
        The routine returns the maximum of the relative error. 
        
        I start at z = 100ish to avoid overflow errors. This test could be 
        constructed better.
        """
        res = True
        try:
            for z in np.linspace(100, 10000.0,10):
                s = sp.sbessel_test_cross_product(z,int(np.floor(z+400)))
                print s
                if s > 1e-13:
                    res = False
        except:
            res = False
        
        self.assertTrue(res)
        
        
        
class TestScalarSphericalTransform(TestCase):
    """Testing the scalar spherical harmonic transforms"""
    def test_with_individual_coefficients_set(self):
        """::Testing spht and ispht with single coefficients, first 10 modes.
        Since these routines haven't been optimized yet, this will take about
        a minute.
        """
        
        res = True
        
        try:
            for n in xrange(0,11):
                for m in xrange(-n,n+1):
                    rnd = np.random.normal(0,10)
                    c = sp.zeros_coefs(48,48)
                    c[n,m] = rnd
                    p = sp.ispht(c,100,100)
                    c2 = sp.spht(p,48,48)
                    s = sp.L2_coef(c - c2)/sp.L2_coef(c)
                    
                    if s > 1e-13:
                        res = False
        except:
            res = False
            
        self.assertTrue(res)
                    
        
