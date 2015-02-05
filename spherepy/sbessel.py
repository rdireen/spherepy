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
             test_sbessel: Spherical Bessel functions 

Randy Direen
1/28/2015


"""
import numpy as np

def sbessely(x,N):
    """Returns a vector of spherical bessel functions yn:

        x:   The argument.
        N:   values of n will run from 0 to N-1.

    """

    out = np.zeros(N,dtype=np.float64)

    out[0] = -np.cos(x)/x
    out[1] = -np.cos(x)/(x**2) - np.sin(x)/x

    for n in xrange(2,N):
        out[n] = ((2.0*n - 1.0) / x)*out[n-1] - out[n-2]

    return out

def sbesselj(x,N):
    """Returns a vector of spherical bessel functions jn:

        x:   The argument.
        N:   values of n will run from 0 to N-1.

    """

    nmax = N-1;
    out = np.zeros(N,dtype=np.float64)
    z = x**2

    out[0] = np.sin(x)/x
    j1 = np.sin(x)/z - np.cos(x)/x

    u = 1
    v = x/(2.0*nmax + 1.0)
    w = v
    n = nmax

    while(np.abs(v/w) > 1e-20):
        n = n+1
        u = 1/(1 - z * u /(4.0 * n**2 - 1.0))
        v *= u - 1
        w += v

    out[nmax] = w

    for n in xrange(nmax-1,0,-1):
        out[n] = 1.0 / ((2.0 * n + 1.0)/x - out[n+1])

    if(np.abs(out[0]) >= np.abs(j1)):
        out[1] *= out[0]
    else:
        out[1] = j1

    for n in xrange(1,nmax):
        out[n+1] *= out[n]

    return out

def sbesselh1(x,N):
    "Spherical Hankel of the first kind"
    
    jn = sbesselj(x,N)
    yn = sbessely(x,N)

    return jn + 1j*yn

def sbesselh2(x,N):
    "Spherical Hankel of the second kind"

    jn = sbesselj(x,N)
    yn = sbessely(x,N)

    return jn - 1j*yn

"""***************************************************************************
******************************************************************************

            The following routines are useful for testing.

******************************************************************************
***************************************************************************"""

def sbesselj_array(xvec,N):
    """Outputs an array where each column is a vector of sbessel values. This
    is useful for plotting a set of Spherical Bessel Functions:

        A = sbessel.sbessely_array(np.linspace(.1,20,100),40)
        for sb in A:
            plot(sb)
        show()
    """

    first_time = True  
    for x in xvec:
        a = sbesselj(x,N)
        if first_time:
            out = np.array([a])
            first_time = False
        else:
            out = np.concatenate([out,[a]], axis=0)
            
    return out.T 

def sbessely_array(xvec,N):
    """Outputs an array where each column is a vector of sbessel values. This
    is useful for plotting a set of Spherical Bessel Functions:

        A = sbessel.sbessely_array(np.linspace(.1,20,100),40)
        for sb in A:
            plot(sb)
        ylim((-.4,.4))
        show()
    """

    first_time = True  
    for x in xvec:
        a = sbessely(x,N)
        if first_time:
            out = np.array([a])
            first_time = False
        else:
            out = np.concatenate([out,[a]], axis=0)
            
    return out.T 

def sbesselj_sum(z,N):
    """Tests the Spherical Bessel function jn using the sum:

        Inf
        sum  (2*n+1) * jn(z)**2 = 1
        n=0


        z:  The argument.
        N:  Large N value that the sum runs too.

    Note that the sum only converges to 1 for large N value (i.e. N >> z).

    The routine returns the relative error of the assumption.
    """

    b = sbesselj(z,N)
    vvv = 2.0*np.array(range(0,N),dtype=np.float64) + 1.0
    sm = np.sum(np.sort(vvv*(b**2)))
    return np.abs((sm - 1.0) / sm) + np.spacing(1)

def sbessel_test_coss_product(z,N):
    """Uses the cross-product relationship to test the routines: 

        j[n+1]y[n] - j[n]y[n+1] = 1 / (z**2)

    where j and y are sbesselj or sbessely vectors for a particular z. Doing 
    this provides a check to see if the bessel functions are being calculated
    correctly.

    The routine returns the maximum of the relative error.
    """

    y1 = sbessely(z,N+1)
    j1 = sbesselj(z,N+1)

    w = y1[:-1]*j1[1:] - j1[:-1]*y1[1:]
    wp = 1.0 / ( (z*np.ones(N,dtype=np.float64)) ** 2)
    return np.max(np.abs(w - wp)/wp)




