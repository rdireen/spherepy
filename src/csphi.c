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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "csphi.h"

/** 
* @file csphi.c
* @brief Low level routines for spherical transforms
* @author Randy Direen
* @version 0.0.1
* @date 2009-09-09
*/

/**
* This is a MACRO that helps in the ynunm routine
* @brief MACRO for ynunm routine
* @see ynunm()
*/
#define B(n,v) (((n)-(v))*((n)+(v)+1.0))

/** 
* Calculate the starting value y(n|nm) for the recursion formula
*that generates the Fourier coefficients of the spherical 
*harmonic Ynm.
* @brief Calculate y(n|nm) 
* 
* @param n The n index in Ynm
* @param m The m index in Ynm
* 
* @return ynnm The starting value y(n|nm) for the recursion formula
* @see ynunm()
*/
double ynnm(int n, int m)
{
	int pm,k;
	double a,ynnm;

	a=1.0/sqrt(4.0*PI);
	pm = abs(m);

	if(n < pm){
		ynnm=0;
	}
	else if (n==0){
		ynnm=a;
	}
	else{
		ynnm=a;
		for(k=1; k<=n; k++){
			ynnm=sqrt((2.0*k+1)/(8.0*k))*ynnm;
		}
		if(n!=pm){
			for(k=n-1; k>=pm; k--){
			ynnm=sqrt(1.0*(n+k+1.0)/(n-k))*ynnm;
			}
		}
	}
	return ynnm;
}

/** 
* Calculates the Fourier coefficients y(nu|n,m). y must have a
* length greater than n.  
* @brief Calculate y(nu|nm) 
* 
* @param en Index n
* @param em Index m
* @param y[] Array of Fourier coefficients for Ynm
* @param len Length of array y
* 
* @see ynnm()
*/
void ynunm(int en,int em,SFLOAT* y,int len)
{
	int k;

	for(k=0;k<len;k++)
		*(y + k)=0;

	if(abs(em) <= en){
		*(y + en) = ynnm(en,em);
		k=en-2;
		if(k >= 0){
			*(y + k) = (B(en,k+1.0)+B(en,k+2.0) - 4.0*em*em)*(*(y + k + 2))/B(en,k);	
			for(k=k-2;k>=0;k-=2){
				(*(y + k))=((B(en,k+1.0)+B(en,k+2.0)-4.0*em*em)*(*(y + k + 2))-B(en,k+3.0)*(*(y + k + 4)))/B(en,k);
			}
		}
	}
}

/** 
 * @brief Finds the next largest number to S, that can be factored into small primes.
 * 
 * @param S The output will be a number greater than or equal to this
 * 
 * @return A number greater than S, that can be factored by 2,3,5,7
 */
int FindQ(int S)
{
	int A;
	int Q = S;

        A = Q;
	while(A != 1)
	{
		if( A % 2 == 0)
			A = A/2;
		else if( A % 3 == 0)
			A = A/3;
		else if( A % 5 == 0)
			A = A/5;
		else if( A % 7 == 0)
			A = A/7;
		else
		{
			A=Q+1;
			Q=A;
		}
	}

	return Q;
}

/** 
 * @brief Calculate s data to be passed to hkm_fc
 * 
 * @param s  a vector of the s data with length Q
 * @param Q
 * @param Nrows number of rows in fdata
 * @param NcoefMax largest spherical coefficient of interest
 */
void SData(SCOMPLEX* s,int Q,int Nrows, int NcoefMax)
{
	int mm,k;
	int nn = NcoefMax + 1;

	memset(s,0,Q*sizeof(SCOMPLEX));

	mm=floor((SFLOAT)Nrows/2.0);

	if((Nrows % 2)==1 ) //data is odd in rows (theta)
	{
		for(k=mm;k<=mm+nn-1;k++)
			if ((k % 2) == 1)
				s[k-mm] = -I/((SCOMPLEX)k);

		for(k=-mm;k<mm;k++)
			if ((abs(k) %2) == 1)
				s[Q+k-mm] = -I/((SCOMPLEX)k);

	}
	else //data is even in rows
	{
		for(k=mm;k<=mm+nn-1;k++)
			if ((k % 2) == 1)
				s[k-mm] = -I/((SCOMPLEX)k);
		for(k=-mm+1;k<mm;k++)
			if ((abs(k) % 2) == 1)
				s[Q+k-mm] = -I/((SCOMPLEX)k);
	}

}
