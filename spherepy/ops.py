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

                ops: Operators on coefficients and patterns 

Randy Direen
2/11/2015

These operations are applied to the theta-phi data and the coefficients. They
are needed for both calculating the scalar and the vector coefficients.

***************************************************************************"""

#---------------------------------------------------------------------Built-ins
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

#---------------------------------------------------------------------3rd Party
import numpy as np

#TODO: Change all xrange instances to range
#and do a 'from six.moves import range' here
from six.moves import xrange

#==============================================================================
# Operators on Fourier Data
#==============================================================================

def fix_even_row_data_fc(fdata):
    """When the number of rows in fdata is even, there is a subtlety that must
    be taken care of if fdata is to satisfy the symmetry required for further
    processing. For an array length of 6,the data is align as [0 1 2 -3 -2 -1] 
    this routine simply sets the row corresponding to the -3 index equal to 
    zero. It is an unfortunate subtlety, but not taking care of this has
    resulted in answers that are not double precision. This operation should
    be applied before any other operators are applied to fdata."""

    L = fdata.shape[0]
    if np.mod(L, 2) == 0:
        fdata[int(L / 2), :] = 0

def pad_rows_fdata(fdata, fdata_extended):
    
    Le = fdata_extended.shape[0]
    L = fdata.shape[0]

    if Le < L:
        raise Exception("Make sure fdata_extended has equal or more rows" + 
                       " than fdata")

    M = int(np.floor(L / 2))

    fdata_extended[:, :] = 0
    fdata_extended[0:M, :] = fdata[0:M, :]
    fdata_extended[-1:-M:-1, :] = fdata[-1:-M:-1, :]

def sin_fc(fdata):
    """Apply sine in the Fourier domain."""

    nrows = fdata.shape[0]
    ncols = fdata.shape[1]

    M = nrows / 2
    fdata[int(M - 1), :] = 0
    fdata[int(M + 1), :] = 0
    
    work1 = np.zeros([nrows, ncols], dtype=np.complex128)
    work2 = np.zeros([nrows, ncols], dtype=np.complex128)

    work1[0, :] = fdata[-1, :]
    work1[1:, :] = fdata[0:-1, :]

    work2[0:-1] = fdata[1:, :]
    work2[-1, :] = fdata[0, :]

    fdata[:, :] = 1.0 / (2 * 1j) * (work1 - work2)

def divsin_fc(fdata):
    """Apply divide by sine in the Fourier domain."""
    
    nrows = fdata.shape[0]
    ncols = fdata.shape[1]

    L = int(nrows / 2) # Assuming nrows is even, which it should be.
    L2 = L - 2 # This is the last index in the recursion for division by sine.
    
    g = np.zeros([nrows, ncols], dtype = np.complex128)
    g[L2, :] = 2 * 1j * fdata[L - 1, :]

    for k in xrange(L2, -L2, -1):
        g[k - 1, :] = 2 * 1j * fdata[k, :] + g[k + 1, :]

    fdata[:, :] = g

def dtheta_fc(fdata):
    """Apply theta derivative in the Fourier domain."""
    
    nrows = fdata.shape[0]
    ncols = fdata.shape[1]
    B = int(nrows / 2) #As always, we assume nrows and ncols are even
    
    a = list(range(0,B))
    ap = list(range(-B,0))
    a.extend(ap)
    
    dtheta = np.zeros([nrows,ncols],np.complex128)
    
    for k in xrange(0,ncols):
        dtheta[:,k] = a
        
    fdata[:,:] = 1j*dtheta * fdata
    

def dphi_fc(fdata):
    """Apply phi derivative in the Fourier domain."""
    
    nrows = fdata.shape[0]
    ncols = fdata.shape[1]
    
    B= int(ncols / 2) #As always, we assume nrows and ncols are even
    
    a = list(range(0,int(B)))
    ap = list(range(-int(B),0))
    a.extend(ap)
    
    dphi = np.zeros([nrows,ncols],np.complex128)
    
    for k in xrange(0,nrows):
        dphi[k,:] = a
        
    fdata[:,:] = 1j*dphi * fdata
    
def sinLdot_fc(tfdata,pfdata):
    """Apply sin of theta times the L operator to the data in the Fourier 
    domain."""
    
    dphi_fc(tfdata)
    
    sin_fc(pfdata)
    dtheta_fc(pfdata)
    
    return 1j*(tfdata - pfdata)

def L_fc(fdata):
    """Apply L in the Fourier domain."""

    fd = np.copy(fdata)

    dphi_fc(fdata)
    divsin_fc(fdata)

    dtheta_fc(fd)

    return (1j * fdata, - 1j * fd)




    
    
    



    