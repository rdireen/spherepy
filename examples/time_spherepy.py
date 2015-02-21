"""***************************************************************************

              time_spherepy: Time the spht and vspht routines. 

Randy Direen
2/18/2015

I'm using this file to see how fast the routines run for various problem 
sizes. I am not using this file to compare to the pure python versions of the 
routines (see benchmark_c_extensions). 

***************************************************************************"""

import spherepy as sp
import numpy as np
import time
import profile

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

Nmax = 500
Nrows = 1024

XX = range(50,550,50)
est = []
eist = []
evt = []
evit = []

print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
print("             Scalar")
print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")

for nm in XX:
    c = sp.random_coefs(nm, nm)
    
    start_time = time.time()
    p = sp.ispht(c, Nrows, Nrows)
    elapsed = time.time() - start_time
    est.append(elapsed)
    print("") 
    print("***Mode Number***: " + str(nm))
    print("")
    print("     Forward Elapsed: " + str(elapsed))

    start_time = time.time()
    c = sp.spht(p, nm, nm)
    elapsed = time.time() - start_time 
    eist.append(elapsed)
    print("     Backward Elapsed: " + str(elapsed))
    print("")


print("")
print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
print("             Vector")
print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
print("")

for nm in XX:
    vc  = sp.random_coefs(nm, nm, coef_type = sp.vector)

    start_time = time.time()
    p = sp.vispht(vc, Nrows, Nrows)
    elapsed = time.time() - start_time 
    evt.append(elapsed)
    
    print("***Mode Number***: " + str(nm))
    print("")
    print("     Forward Elapsed: " + str(elapsed))

    start_time = time.time()
    c = sp.vspht(p, nm, nm)
    elapsed = time.time() - start_time 
    evit.append(elapsed)
    print("     Backward Elapsed: " + str(elapsed))
    print("")


l1, = plt.plot(np.array(XX), np.array(est),label='spht')
l2, = plt.plot(np.array(XX), np.array(eist),label='ispht')
l3, = plt.plot(np.array(XX), np.array(evt),label='vspht')
l4, = plt.plot(np.array(XX), np.array(evit),label='vispht')
plt.legend(handles=[l1, l2, l3, l4],loc=2)

plt.xlabel('Nmax')
plt.ylabel('E-Time [sec]')
plt.show()



