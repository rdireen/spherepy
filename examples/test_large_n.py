from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import numpy as np
import spherepy.csphi as csphi
import spherepy.pysphi as pysphi
from matplotlib.pylab import *

from six.moves import range  

TEST_YNNM_HDR = False
TEST_YNUNM_HDR = False

if TEST_YNNM_HDR:
    #This will test to see that all the starting values within double precision work
    print("Testing modes within double precision")
    for n in range(0, 1010):
        for m in range(-n, n + 1):
            val1 = csphi.ynnm(n, m)
            (val2, ee) = pysphi.ynnm_hdr(n,m)
            valf = (val2 * (10 ** -ee))
            diff =  (val1 - valf) / val1
            assert diff < 1e-12
        print("n = %d ::ynnm = %.18e  ynnm_hdr = (%.18e , %d) : %.18e" % (n,val1, val2, ee, diff))
    print("SUCCESS")
    
   
#plot(list(range(strt, stp)), np.log(nvec))
#show()

N = 1000
(norm, e1) = pysphi.ynnm_hdr(N,N)
(valn, EE) = pysphi.ynunm_hdr(N,N,N+1)
ex10 = 10 ** (np.array(e1 + EE,dtype=np.double))

if TEST_YNUNM_HDR:
    #This will test to see that all the starting values within double precision work
    print("Testing modes within double precision")
    for n in range(0, 1000):
        for m in range(-n, n + 1):
            (norm, e1) = pysphi.ynnm_hdr(n, m)
            (valn, EE) = pysphi.ynunm_hdr(n, m,n+1)
            ex10 = 10 ** (np.array(e1 + EE,dtype=np.double))

            yn_hdr = valn * norm / ex10
            yn = pysphi.ynunm(n, m, n + 1)

            diff =  np.max(yn - yn_hdr) / np.max(yn)
            assert diff < 1e-12
        print("n = %d :: diff : %.18e" % (n, diff))
    print("SUCCESS")

plot( valn * norm / ex10)
show()