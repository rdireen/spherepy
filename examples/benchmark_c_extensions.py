"""***************************************************************************

         benchmark_c_extentions: Compare c extentions to pure python. 

Randy Direen
2/18/2015

This file shows how much faster the c extensions perform over the pure 
python version. 

***************************************************************************"""

import spherepy as sp
import profile
import time

#TODO: Change all xrange instances to range
#and do a 'from six.moves import range' here
from six.moves import xrange

Nmax = 200
Mmax = Nmax
Nrows = 2000
c = sp.random_coefs(Nmax, Nmax)

N = Nmax + 1;
NModes = N + Mmax * (2 * N - Mmax - 1);

head = """********************************************************************

                Benchmark test to compare python version 
                       of the spherical code to the 
                     c version of the spherical code

********************************************************************"""

data_sec = """    A random set of scalar coefficients has been generated with:

                    Nmax = {0}     and     Mmax = {1}
            
    for a total mode count of 
    
                          Total = {2}  Modes                        
                          """.format(str(Nmax),
                                     str(Mmax),
                                     str(NModes)) 


print(head)
print(data_sec)

#transforms without using the c extensions
print("                  -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
print("                  -=-                               -=-")
print("                  -=-   Python version of the code  -=-")
print("                  -=-                               -=-")
print("                  -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
print("")
print("                  Running ispht_slow(scoef,{0},{1})".format(str(Nrows),
                                                            str(Nrows)))
print("")

start_time = time.time()
p = sp.ispht_slow(c, Nrows, Nrows)
elapsed = time.time() - start_time 
f1 = elapsed

print("                  Elapsed Time: {0}".format(str(elapsed)))
print("")
print("                  -------------------------------------")
print("")
print("                  Running spht_slow(patt,{0},{1})".format(str(Nmax),
                                                      str(Nmax)))

start_time = time.time()
c2 = sp.spht_slow(p, Nmax, Nmax)
elapsed = time.time() - start_time
f2 = elapsed
print("")
print("                  Elapsed Time: {0}".format(str(elapsed)))
print("")
err = sp.L2_coef(c - c2) / sp.L2_coef(c)
print("                  Relative Error:  {0}".format(err))


#transforms using the c extensions
print("")
print("                  -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
print("                  -=-                               -=-")
print("                  -=-     C version of the code     -=-")
print("                  -=-                               -=-")
print("                  -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
print("")
print("                  Running ispht(scoef,{0},{1})".format(str(Nrows),
                                                            str(Nrows)))
print("")

start_time = time.time()
p = sp.ispht(c, Nrows, Nrows)
elapsed = time.time() - start_time 
f1 = f1/elapsed

print("                  Elapsed Time: {0}".format(str(elapsed)))
print("")
print("                  -------------------------------------")
print("")
print("                  Running spht(patt,{0},{1})".format(str(Nmax),
                                                      str(Nmax)))

start_time = time.time()
c2 = sp.spht(p, Nmax, Nmax)
elapsed = time.time() - start_time
f2 = f2/elapsed
print("")
print("                  Elapsed Time: {0}".format(str(elapsed)))
print("")
err = sp.L2_coef(c - c2) / sp.L2_coef(c)
print("                  Relative Error:  {0}".format(err))
print("")
print("")
print("          -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
print("          -=-                                               -=-")
print("               ispht() speed up by factor:    {0}  ".format(str(f1)))
print("          -=-                                               -=-")
print("                spht() speed up by factor:    {0}   ".format(str(f2)))
print("          -=-                                               -=-")
print("          -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")




