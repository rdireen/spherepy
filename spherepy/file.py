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

                file: Loading and saving data files 

Randy Direen
2/18/2015

Some file handling routines.

***************************************************************************"""

import numpy as np
import spherepy as sp

def save_patt(patt, filename):

    nrows = patt.nrows
    ncols = patt.ncols
    frmstr = "{0},{1},{2:.16e},{3:.16e}\n"

    ar = patt.array

    with open(filename, 'w') as f: 
        f.write("{0},{1}\n".format(nrows,ncols))
        for nr in xrange(0,nrows):
            for nc in xrange(0,ncols):
                f.write(frmstr.format(nr,nc,ar[nr,nc].real,
                                            ar[nr,nc].imag))
    

def save_coef(scoef, filename):
    
    nmax = scoef.nmax
    mmax = scoef.mmax

    frmstr = "{0:.16e},{1:.16e}\n"

    L = (nmax+1) + mmax*(2*nmax-mmax+1);

    with open(filename, 'w') as f: 
        f.write("{0},{1}\n".format(nmax,mmax))
        for n in xrange(0,L):
            f.write(frmstr.format(scoef._vec[n].real,
                                  scoef._vec[n].imag))


def load_patt(filename):

    with open(filename) as f: 
        lines = f.readlines()

        lst = lines[0].split(',')

        patt = np.zeros([int(lst[0]),int(lst[1])],
                        dtype = np.complex128)

        lines.pop(0)

        for line in lines:
            lst = line.split(',')
            n = int(lst[0])
            m = int(lst[1])
            re = float(lst[2])
            im = float(lst[3])
            patt[n,m] = re + 1j * im

    return sp.ScalarPatternUniform(patt, doublesphere = False)

def load_fdata(filename):

    with open(filename) as f: 
        lines = f.readlines()

        lst = lines[0].split(',')

        patt = np.zeros([int(lst[0]),int(lst[1])],
                        dtype = np.complex128)

        lines.pop(0)

        for line in lines:
            lst = line.split(',')
            n = int(lst[0])
            m = int(lst[1])
            re = float(lst[2])
            im = float(lst[3])
            patt[n,m] = re + 1j * im

    return np.array(patt, dtype = np.complex128)

def load_coef(filename):
    
    with open(filename) as f: 
        lines = f.readlines()

        lst = lines[0].split(',')

        nmax = int(lst[0])
        mmax = int(lst[1])

        L = (nmax+1) + mmax*(2*nmax-mmax+1);

        vec = np.zeros(L,dtype = np.complex128)
   
        lines.pop(0)

        for n,line in enumerate(lines):
            lst = line.split(',')
            re = float(lst[0])
            im = float(lst[1])
            vec[n] = re + 1j * im

    return sp.ScalarCoefs(vec, nmax, mmax)



