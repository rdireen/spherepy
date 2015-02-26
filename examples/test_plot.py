from __future__ import division
import sys
sys.path.append('../spherepy')
import spherepy as sp


#TODO: Change all xrange instances to range
#and do a 'from six.moves import range' here
from six.moves import xrange



fs = """

c[n, m]
=======

2: {7}  {4}  {2}  {6}  {8} 
1:                {3}  {1}  {5} 
0:                               {0}   
n  -------------  -------------  -------------  -------------  -------------  
       m = -2         m = -1         m = 0          m = 1          m = 2        
"""

def _tiny_rep(c):
    sr = "{0:.2}".format(c)
    if sr[0] == '(':
        sr =  sr[1:-1]
    return sr

c = sp.random_coefs(4,4)
c[0,0] = 1


sa = []
cfit = c[0:2,:]
cvec = cfit._vec

for val in cvec:
    sa.append(_tiny_rep(val))

while len(sa) < 9:
    sa.append("")

for n in range(0,9):
    sa[n] = sa[n].center(13)


print fs.format(sa[0],sa[1],sa[2],sa[3],sa[4],sa[5],sa[6],sa[7],sa[8])

    



c = sp.random_coefs(4,1)
cc = c[0:2,:]
p = sp.ispht(cc,50,60)

#sp.plot_sphere_mag(p)