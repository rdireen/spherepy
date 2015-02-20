"""***************************************************************************

Randy Direen
2/19/2015

This script attempts to trip as many exceptions as it can to see what the 
output of the exceptions are. The long strings below show all the nasty 
things that will trip exceptions. When you run this script you will see 
a list responses from all of the exceptions.

***************************************************************************"""


import spherepy as sp
import numpy as np

print("")

template = """
**Test**: {0}
Executed: {1}
Type:     {2}
Message:  {3}
""" 
test_scalar_coefs = True
test_vector_coefs = True
test_scalar_patt_unif = True
test_vector_patt_unif = True

def pt(n, ex, e):
    print(template.format(str(n), ex, str(type(e)), str(e)))
    
def test_executes(executes):
    for n, ex in enumerate(executes):
        try:
            exec(ex)
            print("")
            print("**Test**: {0}".format(str(n)))
            print("PASSED:   {0}".format(ex))
        except Exception as e:
            pt(n, ex, e)
 
if test_scalar_coefs:
    print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")
    print("")
    print("             ScalarCoefs structure with exceptions")
    print("")
    print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")   
    executes_scoefs = ["""zz = sp.zeros_coefs(11, 13)""",
                """zz = sp.ones_coefs(11, 15)""",
                """zz = sp.random_coefs(11, 13)""",
                """scoef[3,2]""",
                """scoef[3,-4]""",
                """scoef[5,7]""",
                """scoef[2,2] = 8""",
                """scoef[2,3] = 8""",
                """scoef[2,:]""",
                """scoef[2,1:3]""",
                """scoef[:,1]""",
                """scoef[1:,3]""",
                """scoef[5,:] = vec""",
                """scoef[:,-1] = vec""",
                """scoef[:,0] = vec""",
                """scoef[4,:] = vec""",
                """scoef[11,0]""",
                """scoef[12,0]""",
                """scoef[11,-7]""",
                """scoef[11,-8]""",
                """scoef[2,1:3] = vec""",
                """scoef[1:,3] = vec""",
                """a = scoef + scoef2""",
                """a = scoef - scoef2""",
                """a = scoef * scoef2""",
                """a = scoef / scoef2""",
                """a = scoef + scoef3""",
                """a = scoef - scoef3""",
                """a = scoef * scoef3""",
                """a = scoef / scoef3""",
                """a = scoef4 + scoef""",
                """a = scoef4 - scoef""",
                """a = scoef4 * scoef""",
                """a = scoef4 / scoef""",
                """a = scoef4 + 2.1""",
                """a = scoef4 - 2.1""",
                """a = scoef4 * 2.1""",
                """a = scoef4 / 2.1""",
                """a = 4.5 + scoef2""",
                """a = 4.5 - scoef2""",
                """a = 4.5 * scoef2""",
                """a = 4.5 / scoef2""",
                """a = scoef / 0"""]
    
    scoef = sp.zeros_coefs(11, 7)
    scoef2 = sp.random_coefs(11, 7)
    scoef3 = sp.random_coefs(11, 8)
    scoef4 = sp.random_coefs(12, 7)
    vec = np.zeros(11, dtype=np.complex128)
    
    test_executes(executes_scoefs)

if test_vector_coefs:
    print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")
    print("")
    print("             VectorCoefs structure with exceptions")
    print("")
    print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")  
    
    
    executes_vcoefs =["""vz = sp.zeros_coefs(11, 13, coef_type=sp.vector)""",
                """vz = sp.ones_coefs(11, 15, coef_type=sp.vector)""",
                """vz = sp.random_coefs(11, 13, coef_type=sp.vector)""",
                """vz = sp.zeros_coefs(11, 10, coef_type=sp.vector)""",
                """vz = sp.ones_coefs(11, 10, coef_type=sp.vector)""",
                """vz = sp.random_coefs(11, 10, coef_type=sp.vector)""",
                """vcoef[3,2]""",
                """vcoef[3,-4]""",
                """vcoef[5,7]""",
                """vcoef[2,2] = 8""",
                """vcoef[2,3] = 8""",
                """vcoef[2,2] = (8,1j)""",
                """vcoef[2,3] = (9,1)""",
                """vcoef[2,:]""",
                """vcoef[2,1:3]""",
                """vcoef[:,1]""",
                """vcoef[1:,3]""",
                """vcoef[5,:] = vec""",
                """vcoef[:,-1] = vec""",
                """vcoef[:,0] = vec""",
                """vcoef[4,:] = vec""",
                """vcoef[11,0]""",
                """vcoef[12,0]""",
                """vcoef[11,-7]""",
                """vcoef[11,-8]""",
                """vcoef[2,1:3] = vec""",
                """vcoef[1:,3] = vec""",
                """a = vcoef + vcoef2""",
                """a = vcoef - vcoef2""",
                """a = vcoef * vcoef2""",
                """a = vcoef / vcoef2""",
                """a = vcoef + vcoef3""",
                """a = vcoef - vcoef3""",
                """a = vcoef * vcoef3""",
                """a = vcoef / vcoef3""",
                """a = vcoef4 + vcoef""",
                """a = vcoef4 - vcoef""",
                """a = vcoef4 * vcoef""",
                """a = vcoef4 / vcoef""",
                """a = vcoef4 + 2.1""",
                """a = vcoef4 - 2.1""",
                """a = vcoef4 * 2.1""",
                """a = vcoef4 / 2.1""",
                """a = 4.5 + vcoef2""",
                """a = 4.5 - vcoef2""",
                """a = 4.5 * vcoef2""",
                """a = 4.5 / vcoef2""",
                """a = vcoef / 0"""]
    
    vcoef = sp.zeros_coefs(11, 7, coef_type=sp.vector)
    vcoef2 = sp.random_coefs(11, 7, coef_type=sp.vector)
    vcoef3 = sp.random_coefs(11, 8, coef_type=sp.vector)
    vcoef4 = sp.random_coefs(12, 7, coef_type=sp.vector)
    
    test_executes(executes_vcoefs)

if test_scalar_patt_unif:
    print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")
    print("")
    print("           ScalarPatternUniform structure with exceptions")
    print("")
    print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")  
    
    
    executes_spatts =[
                """spt = sp.zeros_patt_uniform(11, 14, patt_type=sp.scalar)""",
                """spt = sp.zeros_patt_uniform(11, 15, patt_type=sp.scalar)""",
                """spt = sp.ones_patt_uniform(12, 12, patt_type=sp.scalar)""",
                """spt = sp.ones_patt_uniform(12, 13, patt_type=sp.scalar)""",
                """spt = sp.random_patt_uniform(11, 10, patt_type=sp.scalar)""",
                """spt = sp.random_patt_uniform(11, 11, patt_type=sp.scalar)""",
                """a = spatt + spatt2""",
                """a = spatt - spatt2""",
                """a = spatt * spatt2""",
                """a = spatt / spatt2""",
                """a = spatt + spatt3""",
                """a = spatt - spatt3""",
                """a = spatt * spatt3""",
                """a = spatt / spatt3""",
                """a = spatt4 + spatt""",
                """a = spatt4 - spatt""",
                """a = spatt4 * spatt""",
                """a = spatt4 / spatt""",
                """a = spatt4 + 2.1""",
                """a = spatt4 - 2.1""",
                """a = spatt4 * 2.1""",
                """a = spatt4 / 2.1""",
                """a = 4.5 + spatt2""",
                """a = 4.5 - spatt2""",
                """a = 4.5 * spatt2""",
                """a = 4.5 / spatt2""",
                """a = spatt / 0"""]
    
    spatt = sp.zeros_patt_uniform(11, 8, patt_type=sp.scalar)
    spatt2 = sp.random_patt_uniform(11, 8, patt_type=sp.scalar)
    spatt3 = sp.random_patt_uniform(11, 6, patt_type=sp.scalar)
    spatt4 = sp.zeros_patt_uniform(12, 6, patt_type=sp.scalar)
    
    
    test_executes(executes_spatts)

if test_vector_patt_unif:
    print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")
    print("")
    print("           VectorPatternUniform structure with exceptions")
    print("")
    print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")  
    
    
    executes_vpatts =[
                """vpt = sp.zeros_patt_uniform(11, 14, patt_type=sp.vector)""",
                """vpt = sp.zeros_patt_uniform(11, 15, patt_type=sp.vector)""",
                """vpt = sp.ones_patt_uniform(12, 12, patt_type=sp.vector)""",
                """vpt = sp.ones_patt_uniform(12, 13, patt_type=sp.vector)""",
                """vpt = sp.random_patt_uniform(11, 10, patt_type=sp.vector)""",
                """vpt = sp.random_patt_uniform(11, 11, patt_type=sp.vector)""",
                """a = vpatt + vpatt2""",
                """a = vpatt - vpatt2""",
                """a = vpatt * vpatt2""",
                """a = vpatt / vpatt2""",
                """a = vpatt + vpatt3""",
                """a = vpatt - vpatt3""",
                """a = vpatt * vpatt3""",
                """a = vpatt / vpatt3""",
                """a = vpatt4 + vpatt""",
                """a = vpatt4 - vpatt""",
                """a = vpatt4 * vpatt""",
                """a = vpatt4 / vpatt""",
                """a = vpatt4 + 2.1""",
                """a = vpatt4 - 2.1""",
                """a = vpatt4 * 2.1""",
                """a = vpatt4 / 2.1""",
                """a = 4.5 + vpatt2""",
                """a = 4.5 - vpatt2""",
                """a = 4.5 * vpatt2""",
                """a = 4.5 / vpatt2""",
                """a = vpatt / 0"""]
    
    vpatt  = sp.zeros_patt_uniform(11, 8, patt_type=sp.vector)
    vpatt2 = sp.random_patt_uniform(11, 8, patt_type=sp.vector)
    vpatt3 = sp.random_patt_uniform(11, 6, patt_type=sp.vector)
    vpatt4 = sp.zeros_patt_uniform(12, 6, patt_type=sp.vector)
    
    
    test_executes(executes_vpatts)







    
    
    
