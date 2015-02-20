# Copyright (C) 2015  Randy Direen <spherepy@direentech.com>
#
# This file is part of SpherePy.
#
# SpherePy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SpherePy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SpherePy.  If not, see <http://www.gnu.org/licenses/>

"""***************************************************************************

                     verify: Verification routines

Randy Direen
2/19/2015

These routines are used for comparing SpherePy to other packages. 
Specifically, these routines were put together to compare SpherePy to the 
code written in Matlab. At some point we should use these routines to 
compare results to the NIST code.

***************************************************************************"""

import spherepy as sp
import numpy as np
import file as fl

def verify_spht(pattfile, scoeffile):

    fpatt = fl.load_patt(pattfile)
    fsc = fl.load_coef(scoeffile)
    sc = sp.spht(fpatt,fsc.nmax,fsc.mmax)

    diff = fsc - sc
    
    return (sp.L2_coef(diff) ,sp.LInf_coef(diff),
            sp.L2_coef(diff) / sp.L2_coef(fsc),
            sp.LInf_coef(diff) / sp.LInf_coef(fsc))

def verify_ispht(pattfile, scoeffile):
    
    fsc = fl.load_coef(scoeffile)
    fpatt = fl.load_patt(pattfile)
    patt = sp.ispht(fsc,fpatt.nrows,fpatt.ncols)

    diff = fpatt - patt

    return (sp.L2_patt(diff),sp.LInf_patt(diff),
            sp.L2_patt(diff) / sp.L2_patt(fpatt),
            sp.LInf_patt(diff) / sp.LInf_patt(fpatt))

def verify_fdata(fdatafile, scoeffile):

    fsc = fl.load_coef(scoeffile)
    fdata = fl.load_fdata(fdatafile)
    ffdata = np.zeros([fdata.shape[0], fdata.shape[1]],dtype = np.complex128)
    sp.csphi.sc_to_fc(ffdata,fsc._vec,fsc.nmax,fsc.mmax)

    diff = np.abs(fdata - ffdata)
    return np.max(diff)

def verify_vspht(pattf1, pattf2, scoeff1, scoeff2):
    
    fpatt1 = fl.load_patt(pattf1)
    fpatt2 = fl.load_patt(pattf2)

    fsc1 = fl.load_coef(scoeff1)
    fsc2 = fl.load_coef(scoeff2)

    vfpatt = sp.VectorPatternUniform(fpatt1,fpatt2)
    vfsc = sp.VectorCoefs(fsc1,fsc2,fsc1.nmax,fsc1.mmax)

    vsc = sp.vspht(vfpatt,fsc.nmax,fsc.mmax)
    
    return sp.L2_coef(vfsc - vsc)

def verify_vispht(pattf1,pattf2, scoeff1, scoeff2):
    
    fpatt1 = fl.load_patt(pattf1)
    fpatt2 = fl.load_patt(pattf2)

    fsc1 = fl.load_coef(scoeff1)
    fsc2 = fl.load_coef(scoeff2)

    vfpatt = sp.VectorPatternUniform(fpatt1,fpatt2)
    vfsc = sp.VectorCoefs(fsc1,fsc2,fsc1.nmax,fsc1.mmax)

    vpatt = sp.vispht(vfsc,fpatt1.nrows,fpatt1.ncols)
    
    return sp.L2_patt(vfpatt - vpatt)
