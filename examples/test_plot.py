import sys
sys.path.append('../spherepy')
import spherepy as sp

c = sp.random_coefs(4,1)
cc = c[0:2,:]
p = sp.ispht(cc,50,60)

sp.plot_sphere_mag(p)