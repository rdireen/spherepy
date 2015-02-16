import sys
sys.path.append('../spherepy')
import spherepy as sp

vpatt = sp.ones_patt_uniform(5, 10,patt_type = sp.vector)

vsc = sp.vspht(vpatt,4,4)

c = vsc[:,0]



