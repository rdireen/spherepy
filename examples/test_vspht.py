import sys
sys.path.append('../spherepy')
import spherepy as sp
import numpy as np

#TODO: Change all xrange instances to range
#and do a 'from six.moves import range' here
from six.moves import xrange

c = sp.random_coefs(2,2)
p = sp.ispht(c,3,6)
c2 = sp.spht(p,2,2)

sp.pretty_coefs(c)
sp.pretty_coefs(c2)

res = True
if (sp.L2_coef(c - c2) / sp.L2_coef(c))  > 1e-13:
    print("Failed")
else:
    print("Win")

a = 1








