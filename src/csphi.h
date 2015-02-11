/* Copyright (C) 2015  Randy Direen <spherepy@direentech.com>
*
* This file is part of SpherePy.
*
* SpherePy is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* SpherePy is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with SpherePy.  If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __CSPHI__
#define __CSPHI__

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "kiss_fft.h"

#define SFLOAT double
#define SINT int

typedef struct {
    SFLOAT r;
    SFLOAT i;
} SCOMPLEX;

#define SUCCESS    0
#define OUTBOUNDS  1

#define PI 3.14159265358979323846



SFLOAT ynnm(int n,int m);

void ynunm(int en,int em,SFLOAT* y,int len);

int FindQ(int S);

void SData(SCOMPLEX* s,int Q,int Nrows, int NcoefMax);

void hkm_fc(SCOMPLEX* gcoef,int Nrow,int Ncol, 
            int n, int m, 
            SCOMPLEX* hkm, int len,
            SCOMPLEX* ss, int Q,
            SCOMPLEX* ff, int Q2,
            kiss_fft_cfg kiss_cfg_fw,
            kiss_fft_cfg kiss_cfg_bw);

void bnm_fc(SCOMPLEX * fdata,int Nrow, int Ncol, 
            int Nmax, int m,
            SCOMPLEX* vec, int L,
            SCOMPLEX* ss,int Q,
            SCOMPLEX* ff, int Q2,
            SCOMPLEX* hkm, int Lhkm,
            SFLOAT* y, int Ly,
            kiss_fft_cfg kiss_cfg_fw,
            kiss_fft_cfg kiss_cfg_bw);

void fc_to_sc(SCOMPLEX* fdata, int Nrow, int Ncol,
              SCOMPLEX* sc, int L,
              int Nmax, int Mmax);

#endif //__CSPHI__
