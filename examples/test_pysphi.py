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


import spherepy as sp
import numpy as np

T1 = False
T2 = True
T3 = False
T4 = True
T5 = True

if T1:
    nrows = 50
    ncols = 50
    out = sp.pysphi.sph_harmonic_tp(nrows,ncols,5,2)
    T = np.abs(out)
    sp.plot_mag_on_sphere(T)

if T2:
    nrows = 10
    Nmax = 2
    Q = sp.pysphi.smallest_prime_factor(nrows + Nmax )
    sd = sp.pysphi.s_data(nrows, Nmax, Q)

if T3:

    a = [[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9]]
    fdata = np.array(a,dtype=np.complex64)

    nrows = 10
    Nmax = 3
    Q = sp.pysphi.smallest_prime_factor(nrows + Nmax - 1)
    sd = sp.pysphi.s_data(nrows, Nmax, Q)

    b = sp.pysphi.hkm_fc(fdata,Nmax,-1,sd)

if T4:

    a = [[1, 2, 3, 4, 5, 6, 7, 8],
         [3, 2, 3, 2, 3, 5, 9, 5],
         [1, 5, 3, 1, 5, 3, 1, 5]]

    b = [[1, 2, 3, 4, 5, 6, 7, 8],
         [3, 2, 3, 2, 3, 5, 9, 5],
         [1, 5, 3, 1, 5, 3, 1, 5],
         [2, 3, 2, 3, 2, 1, 1, 1]]

    ar1 = np.array(a,dtype=np.complex128)
    ar2 = np.array(b,dtype=np.complex128)
    d1 = sp.continue_sphere(ar1,1.0)
    d2 = sp.continue_sphere(ar2,1.0)

    ssphere1 = sp.ScalarPatternUniform(d1)
    ssphere2 = sp.ScalarPatternUniform(ar2)

    print(ssphere1.single_val)
    print(ssphere2.single_val)
    print(ssphere1.is_symmetric)
    print(ssphere2.is_symmetric)


    vvv1 = sp.pysphi.fc_to_sc(ssphere2._dsphere,2,2)

    #pysphi.fix_even_row_data_fc(d2)

    #fdata_extended = np.zeros([8,8],dtype=np.complex128)

    #pysphi.pad_rows_fdata(d2,fdata_extended)

    #pysphi.sin_fc(fdata_extended)

    #this matches Matlab
    scoef = sp.spht(ssphere2,2,2)

    th = scoef._array_2d_repr()

    rs = scoef._reshape_n_vecs()

    #ps.pcolor_coefs(scoef)
    sp.plot_coefs(scoef)

    print scoef[:,2]

    sp = sp.ispht(scoef,6,8)

if T5:
    rr = sp.random_coefs(10,10) + 1j*sp.random_coefs(10,10)

    ss = 1j - 2 / (sp.zeros_coefs(5,4) + .001) + \
            4* sp.random_coefs(5,4) / 6.0
    ss += sp.ones_coefs(5,4,coef_type=sp.scalar)

    qq = 1 + 2*sp.ones_coefs(3,3) / 4*sp.ones_coefs(3,3) - 2

    ds1 =1j + 4*sp.ones_patt_uniform(5,8)/3 - \
        sp.ones_patt_uniform(5,8) -1 + \
        10*sp.random_patt_uniform(5,8)

    print ds1.single_val

    print sp.array(ds1)

    zc = sp.zeros_coefs(10,10)
    zc[2,1] = 1.0
    zc[2,-2] = 3.0
    zc[0,0] = 1j*5.0
    zc[5,1] = 7.0

    print zc[2,1]
    print zc[2,-2]
    print zc[5,1]

    




vv = raw_input()







