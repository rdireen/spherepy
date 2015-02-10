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
#include <complex.h>
#include <string.h>

#define SFLOAT double
#define SCOMPLEX double complex
#define SINT int

#define SUCCESS    0
#define OUTBOUNDS  1

#define PI 3.14159265358979323846

SFLOAT ynnm(int n,int m);
void ynunm(int en,int em,SFLOAT* y,int len);
int FindQ(int S);
void SData(SCOMPLEX* s,int Q,int Nrows, int NcoefMax);

#endif //__CSPHI__
