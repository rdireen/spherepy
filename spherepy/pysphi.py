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


"""
             pyphi: Low Level Routines 

Randy Direen
1/28/2015


"""
import numpy as np

def ynnm(n, m):
    """Initial value for recursion formula""" 
    a = 1.0 / np.sqrt(4.0 * np.pi)
    pm = np.abs(m)

    out = 0.0

    if(n < pm):
        out = 0.0
    elif(n == 0):
        out = a
    else:
        out = a
        for k in xrange(1, n + 1):
            out *= np.sqrt((2.0 * k + 1.0) / 8.0 / k)

        if(n != pm):
            for k in xrange(n - 1, pm - 1, -1):
                out *= np.sqrt((n + k + 1.0) / (n - k))
    return out

def ynunm(n, m, L):
    """Fourier coefficients for spherical harmonics"""

    out = np.zeros(L, dtype=np.float64)
    tmp1 = 0 
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0       
    if(np.abs(m) <= n):
        out[n] = ynnm(n, m)
        k = n - 2
        if(k >= 0):
            tmp1 = (n - k - 1.0) * (n + k + 2.0)
            tmp2 = (n - k - 2.0) * (n + k + 3.0) - 4.0 * m ** 2
            tmp4 = ((n - k) * (n + k + 1.0))
            out[k] = (tmp1 + tmp2) * out[k + 2] / tmp4

            for k in xrange(n - 4, -1, -2):
                tmp1 = (n - k - 1.0) * (n + k + 2.0)
                tmp2 = (n - k - 2.0) * (n + k + 3.0) - 4.0 * m ** 2
                tmp3 = (n - k - 3.0) * (n + k + 4.0);
                tmp4 = ((n - k) * (n + k + 1.0))
                out[k] = ((tmp1 + tmp2) * out[k + 2] - tmp3 * out[k + 4]) / tmp4
    return out    

def ynunm_work(n, m, work):
    """Fourier coefficients for spherical harmonics"""

    work[:] = 0
    tmp1 = 0 
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0       
    if(np.abs(m) <= n):
        work[n] = ynnm(n, m)
        k = n - 2
        if(k >= 0):
            tmp1 = (n - k - 1.0) * (n + k + 2.0)
            tmp2 = (n - k - 2.0) * (n + k + 3.0) - 4.0 * m ** 2
            tmp4 = ((n - k) * (n + k + 1.0))
            work[k] = (tmp1 + tmp2) * work[k + 2] / tmp4

            for k in xrange(n - 4, -1, -2):
                tmp1 = (n - k - 1.0) * (n + k + 2.0)
                tmp2 = (n - k - 2.0) * (n + k + 3.0) - 4.0 * m ** 2
                tmp3 = (n - k - 3.0) * (n + k + 4.0);
                tmp4 = ((n - k) * (n + k + 1.0))
                work[k] = ((tmp1 + tmp2) * work[k + 2] - tmp3 * work[k + 4]) / tmp4
    return work
    
def sph_harmonic_tp(nrows, ncols, n, m):
    """Produces Ynm(theta,phi)

        theta runs from 0 to pi and has 'nrows' points.

        phi runs from 0 to 2*pi - 2*pi/ncols and has 'ncols' points.
    
    """
    nuvec = ynunm(n, m, n + 1)
    num = nuvec[1:] * (-1) ** m
    num = num[::-1]
    yvec = np.concatenate([num, nuvec])
    
    theta = np.linspace(0.0, np.pi, nrows)
    phi = np.linspace(0.0, 2.0 * np.pi - 2.0 * np.pi / ncols)
    nu = np.array(range(-n, n + 1), dtype=np.complex128)
    
    out = np.zeros([nrows, ncols], dtype=np.complex128)

    for nn in xrange(0, nrows):
        exp_vec = np.exp(1j * nu * theta[nn])
        for mm in xrange(0, ncols):          
            vl = np.sum(yvec * exp_vec)
            out[nn, mm] = 1j ** m * np.exp(1j * m * phi[mm]) * vl

    return out

def smallest_prime_factor(Q):
    """Find the smallest number, factorable by the small primes 2,3,4, and 7 
    that is larger than the argument Q"""

    A = Q;
    while(A != 1):
        if(np.mod(A, 2) == 0):
            A = A / 2
        elif(np.mod(A, 3) == 0):
            A = A / 3
        elif(np.mod(A, 5) == 0):
            A = A / 5
        elif(np.mod(A, 7) == 0):
            A = A / 7;
        else:
            A = Q + 1;
            Q = A;

    return Q

def s_data(nrows_fdata, Nmax, Q):
    """ I am going to assume we will always have even data. This is pretty 
    safe because it means that we have measured both poles of the sphere and 
    have data that has been continued.

        nrows_fdata:  Number of rows in fdata.
        Nmax:         The largest number of n values desired.
        Q:            A value greater than nrows_fdata + Nmax. This can be
                      selected to be factorable into small primes to 
                      increase the speed of the fft (probably not that big 
                      of a deal today).

    """

    if np.mod(nrows_fdata, 2) == 1:
        raise Exception("nrows_fdata must be even.")
    
    L1 = nrows_fdata

    s = np.zeros(Q, dtype=np.complex128)
    MM = L1 / 2

    for nu in xrange(-MM, MM + Nmax + 1):
        if np.mod(nu, 2) == 1:
            s[nu - MM] = -1j / nu

    return s

def hkm_fc(fdata, Nmax, m, s):
    """ Assume fdata has even rows"""

    f = fdata[:, m]
    L1 = f.size
    MM = L1 / 2
    Q = s.size

    ff = np.zeros(Q, dtype=np.complex128)
    for n in xrange(MM, L1):
        ff[n] = f[n - MM]

    for n in xrange(0, MM):
        ff[n] = f[n + MM]

    # For larger problems, this speeds things up pretty good.
    F = np.fft.fft(ff)
    S = np.fft.fft(s)
    out = 4 * np.pi * np.fft.ifft(F * S)

    return out[0:Nmax + 1]

def bnm_vec_fc(fdata, Nmax, m):

    nrows = fdata.shape[0]
    Q = smallest_prime_factor(nrows + Nmax)
    s = s_data(nrows, Nmax, Q)
    h = hkm_fc(fdata, Nmax, m, s)

    absm = np.abs(m)

    out = np.zeros(Nmax - absm + 1, dtype=np.complex128)

    for n in xrange(absm, Nmax + 1):

        ynm = ynunm(n, m, n + 1)

        out[n - absm] = 1j ** (-m) * h[0] * ynm[0]

        if n != 0:
            out[n - absm] += 1j ** (-m) * 2 * np.sum(h[1:n + 1] * ynm[1:n + 1])
 
    return out


def bnm_vec_fc_work(fdata, Nmax, m, work):

    nrows = fdata.shape[0]
    Q = smallest_prime_factor(nrows + Nmax)
    s = s_data(nrows, Nmax, Q)
    h = hkm_fc(fdata, Nmax, m, s)

    absm = np.abs(m)

    out = np.zeros(Nmax - absm + 1, dtype=np.complex128)

    for n in xrange(absm, Nmax + 1):

        ynunm_work(n, m, work)

        out[n - absm] = 1j ** (-m) * h[0] * work[0]

        if n != 0:
            out[n - absm] += 1j ** (-m) * 2 * np.sum(h[1:n + 1] * work[1:n + 1])
 
    return out

def fc_to_sc(gcoef, Nmax, Mmax):
    
    c = bnm_vec_fc(gcoef, Nmax, 0)
    
    for m in xrange(1, Mmax + 1):
        a = bnm_vec_fc(gcoef, Nmax, -m)  
        c = np.concatenate([c, a]) 
        a = bnm_vec_fc(gcoef, Nmax, m)  
        c = np.concatenate([c, a])

    return c

def fc_to_sc_work(gcoef, Nmax, Mmax, work):
    
    c = bnm_vec_fc_work(gcoef, Nmax, 0, work)
    
    for m in xrange(1, Mmax + 1):
        a = bnm_vec_fc_work(gcoef, Nmax, -m, work)  
        c = np.concatenate([c, a]) 
        a = bnm_vec_fc_work(gcoef, Nmax, m, work)  
        c = np.concatenate([c, a])

    return c

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
        fdata[L / 2, :] = 0

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
    """Apply sin in the Fourier domain"""

    nrows = fdata.shape[0]
    ncols = fdata.shape[1]

    M = nrows / 2
    fdata[M - 1, :] = 0
    fdata[M + 1, :] = 0
    
    work1 = np.zeros([nrows, ncols], dtype=np.complex128)
    work2 = np.zeros([nrows, ncols], dtype=np.complex128)

    work1[0, :] = fdata[-1, :]
    work1[1:, :] = fdata[0:-1, :]

    work2[0:-1] = fdata[1:, :]
    work2[-1, :] = fdata[0, :]

    fdata[:, :] = 1.0 / (2 * 1j) * (work1 - work2)

def mindx(m, nmax, mmax):
    """index to the first n value for a give m within the spherical 
    coefficients vector. Used by sc_to_fc"""

    ind = 0
    NN = nmax + 1

    if np.abs(m) > mmax:
        raise Exception("|m| cannot be larger than mmax")

    if (m != 0):
        ind = NN
        ii = 1
        for i in xrange(1, np.abs(m)):
            ind = ind + 2 * (NN - i)
            ii = i + 1

        if m > 0:
            ind = ind + NN - ii

    return ind

def fcvec_m_sc(vec, m, nmax, nrows):
    
    F = np.zeros(nrows, dtype=np.complex128)
    K = nmax + 1 

    for n in xrange(np.abs(m), K):
        ynm = ynunm(n, m, K)
        F[0:nmax + 1] += vec[n - np.abs(m)] * ynm.T

    F[0:K] = F[0:K] * 1j ** m

    mm = (-1) ** m
    if np.mod(nrows, 2) == 0:
        H = nrows / 2 - 1
    else:
        H = (nrows - 1) / 2

    for k in xrange(0, H):
        F[-(k + 1)] = mm * F[k + 1]

    return F

def sc_to_fc(spvec, nmax, mmax, nrows, ncols):
    """assume Ncols is even"""

    if np.mod(ncols, 2) == 1:
        raise Exception("ncols is required to be even")

    fdata = np.zeros([nrows, ncols], dtype=np.complex128)

    for k in xrange(0, ncols / 2):
        if k < mmax:
            kk = k
            ind = mindx(kk, nmax, mmax)
            vec = spvec[ind:ind + nmax - np.abs(kk) + 1]
            fdata[:, kk] = fcvec_m_sc(vec, kk, nmax, nrows)

            kk = -(k + 1)
            ind = mindx(kk, nmax, mmax)
            vec = spvec[ind:ind + nmax - np.abs(kk) + 1]
            fdata[:, kk] = fcvec_m_sc(vec, kk, nmax, nrows)

        if k == mmax:
            kk = k
            ind = mindx(kk, nmax, mmax)
            vec = spvec[ind:ind + nmax - np.abs(kk) + 1]
            fdata[:, kk] = fcvec_m_sc(vec, kk, nmax, nrows)

    return fdata   


def fcvec_m_sc_work(vec, m, nmax, nrows, work):
    
    F = np.zeros(nrows, dtype=np.complex128)
    K = nmax + 1 

    for n in xrange(np.abs(m), K):
        ynunm_work(n, m, work)
        F[0:nmax + 1] += vec[n - np.abs(m)] * work[0:nmax + 1].T

    F[0:K] = F[0:K] * 1j ** m

    mm = (-1) ** m
    if np.mod(nrows, 2) == 0:
        H = nrows / 2 - 1
    else:
        H = (nrows - 1) / 2

    for k in xrange(0, H):
        F[-(k + 1)] = mm * F[k + 1]

    return F

def sc_to_fc_work(spvec, nmax, mmax, nrows, ncols, work):
    """assume Ncols is even"""

    if np.mod(ncols, 2) == 1:
        raise Exception("ncols is required to be even")

    fdata = np.zeros([nrows, ncols], dtype=np.complex128)

    for k in xrange(0, ncols / 2):
        if k < mmax:
            kk = k
            ind = mindx(kk, nmax, mmax)
            vec = spvec[ind:ind + nmax - np.abs(kk) + 1]
            fdata[:, kk] = fcvec_m_sc_work(vec, kk, nmax, nrows, work)

            kk = -(k + 1)
            ind = mindx(kk, nmax, mmax)
            vec = spvec[ind:ind + nmax - np.abs(kk) + 1]
            fdata[:, kk] = fcvec_m_sc_work(vec, kk, nmax, nrows, work)

        if k == mmax:
            kk = k
            ind = mindx(kk, nmax, mmax)
            vec = spvec[ind:ind + nmax - np.abs(kk) + 1]
            fdata[:, kk] = fcvec_m_sc_work(vec, kk, nmax, nrows, work)

    return fdata   


   
