"""
Randy Direen
2/21/2015

Profiling script to see where the bottlenecks are in the spht, vspht, ispht, 
and the vispht routines.  Set ctrl below, to one of the tests.

"""
import spherepy as sp
import numpy as np
import profile

#TODO: Change all xrange instances to range
#and do a 'from six.moves import range' here
from six.moves import xrange

Nmax = 400
Nrows = 1024

ISPHT = 0
SPHT = 1
VISPHT = 2
VSPHT = 3
FFT = 4

ctrl = SPHT

c = sp.random_coefs(Nmax, Nmax)
vc  = sp.random_coefs(Nmax, Nmax, coef_type = sp.vector)

if ctrl == ISPHT:
    profile.run('p = sp.ispht(c,Nrows,Nrows)',sort=1)

if ctrl == SPHT:
    p = sp.ispht(c, Nrows, Nrows)
    profile.run('c2 = sp.spht(p,Nmax,Nmax)',sort=1)

if ctrl == VISPHT:
    profile.run('p = sp.vispht(vc,Nrows,Nrows)',sort=1)

if ctrl == VSPHT:
    vp = sp.vispht(vc, Nrows, Nrows)
    profile.run('vc2 = sp.vspht(vp,Nmax,Nmax)',sort=1)

if ctrl == FFT:
    """Interesting, fft2 calls cfft twice"""

    data = np.ones([Nrows,Nrows],dtype = np.complex128)

    profile.run('fdata = np.fft.fft2(data) / (Nrows * Nrows)',sort=1)



