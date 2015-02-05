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
             test_sbessel: Spherical Bessel functions test

Randy Direen
1/28/2015


"""

import pysphi
import sbessel
import plot_sphere as ps
import numpy as np
import random as rn
from pylab import *
print("Test Spherical Bessel Routines")

T1 = False
T2 = False
T3 = False
T4 = False

#*****************************************************************************
#Test the jn sbessel function summation
if T1:
    N = 10000;
    jsum = np.zeros(N,dtype=np.float64)
    xvec = np.zeros(N,dtype=np.float64)

    for n,x in enumerate(np.linspace(0.001, N,400)):
        xx =  x + rn.uniform(-10,10)
        xvec[n] = xx
        jsum[n] = sbessel.sbesselj_sum(xx,int(x)+400)

    figure(1)
    with np.errstate(divide='ignore'):
        plot(xvec,10*np.log10(jsum))

#*****************************************************************************
#Test the cross-product relationship
if T2:
    N = 100
    bbvec = np.zeros(N,dtype=np.float64)
    for n,x in enumerate(np.linspace(1,50,N)):
        bbvec[n] = sbessel.sbessel_test_coss_product(x,100)

    figure(2)
    plot(10*np.log10(bbvec))


#*****************************************************************************
if T3:
    figure(3)
    x = np.linspace(.1,20,100)
    A = sbessel.sbesselj_array(x,5)
    for bs in A:
        plot(x,bs)
    ylim((-1,1))
    xlabel("x")
    ylabel("jn")
    title("First few sbesselj functions")

#*****************************************************************************
if T4:
    figure(4)
    x = np.linspace(.1,20,100)
    A = sbessel.sbessely_array(x,5)
    for bs in A:
        plot(x,bs)
    ylim((-.4,.4))
    xlabel("x")
    ylabel("yn")
    title("First few sbesselj functions")

show()
