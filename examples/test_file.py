import sys
sys.path.append('../spherepy')
import spherepy as sp

#TODO: Change all xrange instances to range
#and do a 'from six.moves import range' here
from six.moves import xrange


patt1 = sp.random_patt_uniform(5, 10,patt_type = sp.scalar)

sp.file.save_patt(patt1,'t1.txt')

patt2 = sp.file.load_patt('t1.txt')

a = sp.L2_patt(patt1 - patt2)

scoef1 = sp.random_coefs(10,7)

sp.file.save_coef(scoef1,'t2.txt')

scoef2 = sp.file.load_coef('t2.txt')

a = sp.L2_coef(scoef1 - scoef2)

a = 1





