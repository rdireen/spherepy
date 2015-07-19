from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import numpy as np
import spherepy.csphi as csphi
import spherepy.pysphi as pysphi
from matplotlib.pylab import *

from six.moves import range  

strt = 1050
stp = 1100

vec = []
for n in range(1050, 1100):
    val = csphi.ynnm(n,5)
    vec.append(val)
    print("n = %d  ynnm(n,0) = %.16f" % (n, val))

nvec = np.array(vec, dtype=np.double)
plot(list(range(strt, stp)), np.log(nvec))
show()