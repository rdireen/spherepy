import sys
sys.path.append('../spherepy')
import spherepy as sp

vpatt = sp.ones_patt_uniform(5, 10,patt_type = sp.vector)

vsc = sp.vspht(vpatt,4,4)

c = vsc[:,0]

c = sp.zeros_coefs(100, 100)
c[2, 0] = 1.0
p = sp.ispht(c, 102, 202)
c2 = sp.spht(p, 100, 100)



