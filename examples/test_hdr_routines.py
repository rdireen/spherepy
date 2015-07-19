from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import numpy as np
import spherepy.csphi as csphi
import spherepy.pysphi as pysphi
from matplotlib.pylab import *
import spherepy as sp

Nrows = 100
Ncols = 202
patt = sp.random_patt_uniform(Nrows,Ncols)
c = sp.spht_ldr(patt, Nrows - 2, Nrows -2)
ch = sp.spht(patt, Nrows - 2, Nrows -2)

print(sp.L2_coef(c -ch))

#sp.pcolor_coefs(c - ch)
#rc = sp.random_coefs(Nrows- 2, Nrows - 2)
p = sp.ispht_ldr(c, Nrows, Ncols)
ph = sp.ispht(c, Nrows, Ncols)

pcolormesh((np.abs(p.array - ph.array)))
print(np.max(np.abs(p.array - ph.array)))
show()
