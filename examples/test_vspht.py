import sys
sys.path.append('../spherepy')
import spherepy as sp
import numpy as np

#TODO: Change all xrange instances to range
#and do a 'from six.moves import range' here
from six.moves import xrange

vpatt = sp.ones_patt_uniform(5, 10,patt_type = sp.vector)

vsc = sp.vspht(vpatt,4,4)

c = vsc[:,0]

c = sp.zeros_coefs(100, 100)
c[2, 0] = 1.0
p = sp.ispht(c, 102, 202)
c2 = sp.spht(p, 100, 100)


a = [[1, 2, 3, 4, 5, 6, 7, 8],
     [3, 2, 3, 2, 3, 5, 9, 5],
     [1, 5, 3, 1, 5, 3, 1, 5],
     [2, 3, 2, 3, 2, 1, 1, 1]]

b = [[1, 5, 3, 4, 5, 5, 7, 8],
     [3, 5, 3, 3, 3, 5, 9, 5],
     [1, 5, 3, 2, 5, 5, 1, 5],
     [2, 5, 2, 1, 2, 5, 1, 1]]

aa = np.array(a,dtype = np.complex128)
bb = np.array(b,dtype = np.complex128)
vv = sp.VectorPatternUniform(aa,bb)

vcf = sp.vspht(vv,3,3)
vpatt = sp.vispht(vcf, 4, 8)

nv = vcf[:,1]

a = 1



