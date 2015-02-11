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
#include <stdbool.h>

#include "csphi.h"
#include "kiss_fft.h"

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
 * @param Nmax largest spherical coefficient of interest
 */
void SData(SCOMPLEX* s,int Q,int Nrows, int Nmax)
{
	int mm,k;
	int nn = Nmax + 1;

	memset(s,0,Q*sizeof(SCOMPLEX));

        //we assume Nrows is even
	mm=(SFLOAT)Nrows/2.0;

	for(k=mm;k<=mm+nn - 1;k++)
		if ((k % 2) == 1)
			s[k-mm].i = -1/((SFLOAT)k);
	for(k=-mm;k<mm;k++)
		if ((abs(k) % 2) == 1)
			s[Q+k-mm].i = -1/((SFLOAT)k);
}

/**
 * Calculates hkm which is an operation on the fourier coefficients (fc).
 * The resulting output hkm is used to calculate the coefficients bnm. The 
 * parameter Q must be greater than Nrow+n and there is no error handling to 
 * correct for Q any smaller.
 * @brief Calculate hkm
 * 
 * @param gcoef Fourier coefficients
 * @param Nrow  Number of rows in gcoef
 * @param Ncol  Number of cols in gcoef
 * @param n     bnm n index (when used int bnm, this is tha maximum n value)
 * @param m     bnm m index
 * @param hkm[]	Output 
 * @param len   Length of output, must be greater than n and less than or equal to Q (set to n+1)
 * @param Q     A number larger than Nrow+n, make factorable into small primes for speed
 */
void hkm_fc(SCOMPLEX* gcoef,int Nrow,int Ncol, 
            int n, int m, 
            SCOMPLEX* hkm, int len,
            SCOMPLEX* ss, int Q,
            SCOMPLEX* ff, int Q2,
            kiss_fft_cfg kiss_cfg_fw,
            kiss_fft_cfg kiss_cfg_bw)
{
	int k, ind, mm;

	memset(ff,0,Q*sizeof(SCOMPLEX));
	memset(hkm,0,len*sizeof(SCOMPLEX));

	// properly get the index into gcoef that m corresponds to
	if(m>=0)
		ind = m;
	else
		ind = Ncol+m;
	
        // we assume the number of rows is even
	mm = Nrow / 2;

	for(k=0;k<mm;k++)
	{
		ff[k] = gcoef[Ncol*(mm+k)+ind];
		ff[mm+k] = gcoef[Ncol*k+ind];
	}

        kiss_fft(kiss_cfg_fw,(kiss_fft_cpx*)ff,(kiss_fft_cpx*)ff);        

	for(k=0;k<Q;k++)
        {
		ff[k].r = 4.0*PI*(ss[k].r*ff[k].r - ss[k].i*ff[k].i);
                ff[k].i = 4.0*PI*(ss[k].r*ff[k].i + ss[k].i*ff[k].r);
        }

	kiss_fft(kiss_cfg_bw,(kiss_fft_cpx*)ff,(kiss_fft_cpx*)ff);

	for(k=0; k<len; k++)
        {
		hkm[k].r = ff[k].r / (SFLOAT)Q;
                hkm[k].i = ff[k].i / (SFLOAT)Q;
        } 

}

/** 
 * @brief Calcutates a vector vec of spherical coefficients
 * 
 * @param fdata Fourier coefficient data
 * @param Nrow  Number of rows in the Fourier coefficient data
 * @param Ncol  Number of columns in the Fourier coefficient data
 * @param Nmax Maximum number of spherical coefficients desired
 * @param m  The m index into the spherical coefficients
 * @param vec The output vector of length Nmax-abs(m)+1
 * @param s The s data given by SData
 * @param Q The length of s
 */
void bnm_fc(SCOMPLEX * fdata,int Nrow, int Ncol, 
            int Nmax, int m,
            SCOMPLEX* vec, int L,
            SCOMPLEX* ss,int Q,
            SCOMPLEX* ff, int Q2,
            SCOMPLEX* hkm, int Lhkm,
            SFLOAT* y, int Ly,
            kiss_fft_cfg kiss_cfg_fw,
            kiss_fft_cfg kiss_cfg_bw)
{
	int n,k;
	int absm = abs(m);

	hkm_fc(fdata, Nrow, Ncol, 
               Nmax, m,
               hkm, Lhkm,
               ss, Q, 
               ff, Q2,
               kiss_cfg_fw,
               kiss_cfg_bw);

	for(n=absm;n<=Nmax;n++)
	{
		ynunm(n,m,y,Ly);
		
		vec[n-absm].r = hkm[0].r*y[0];
                vec[n-absm].i = hkm[0].i*y[0];

		for(k = 1;k<n+1;k++)
		{
                        //TODO: need to sort before I sum
			vec[n-absm].r += 2.0*hkm[k].r*y[k];
                        vec[n-absm].i += 2.0*hkm[k].i*y[k];
		}

                //In the following: vec[n-absm] = i**(-m) * vec[n-absm] 
                if(m % 2 == 1)
                {
		        vec[n-absm].r = pow(-1, (m-1)/2) * vec[n-absm].i;
                        vec[n-absm].i = -pow(-1, (m-1)/2) * vec[n-absm].r;
                }
                else
                {
                        vec[n-absm].r = pow(-1, m/2) * vec[n-absm].r;
                        vec[n-absm].i = pow(-1, m/2) * vec[n-absm].i;

                }
	}
}


/** 
 * @brief Fourier coefficients to spherical harmonic coefficients
 * 
 * @param fdata Fourier coefficients
 * @param Nrow Number of rows in fdata
 * @param Ncol Number of columns in fdata
 * @param sc Array of coefficients
 */
void fc_to_sc(SCOMPLEX* fdata, int Nrow, int Ncol,
              SCOMPLEX* sc, int L,
              int Nmax, int Mmax)
{

	int m, Q;	
        SCOMPLEX* pt;
        int* inds;
        int N = Nmax + 1;
        int QQ = Nmax;

	Q = FindQ(Nrow+Nmax);

        //Allocate all memory here
	SCOMPLEX* s = (SCOMPLEX *)malloc(Q*sizeof(SCOMPLEX));
        SData(s,Q, Nrow,Nmax);

        kiss_fft_cfg kiss_cfg_fw = kiss_fft_alloc(Q,false,0,0);
        kiss_fft_cfg kiss_cfg_bw = kiss_fft_alloc(Q,true,0,0);

	kiss_fft(kiss_cfg_fw,(kiss_fft_cpx*)s,(kiss_fft_cpx*)s);

        SCOMPLEX* hkm = (SCOMPLEX*)malloc((Nmax+1)*sizeof(SCOMPLEX));
	SFLOAT* y = (SFLOAT*)malloc((Nmax+1)*sizeof(SFLOAT));	
        SCOMPLEX* ff = (SCOMPLEX*)malloc(Q*sizeof(SCOMPLEX));

        inds = (int*)malloc((2*Mmax+1)*sizeof(int));
	memset(inds,0,(2*Mmax+1)*sizeof(int));


        //setup indices that index into sc
        inds[0] = 0;
	for(m=1; m<=Mmax; m++)
	{
		inds[2*m-1] = N;
		N = N + QQ;
		inds[2*m] = N;
		N = N + QQ;
		QQ--;
	}


        //Calculate each of the bnm coefficients
	bnm_fc(fdata, Nrow, Ncol,
               Nmax, 0,
               sc, Nmax + 1,
               s, Q,
               ff, Q,
               hkm, Nmax+1,
               y, Nmax+1,
               kiss_cfg_fw,
               kiss_cfg_bw);

	for(m=1; m <= Mmax; m++)
	{
		pt = sc + inds[2*m-1];
		bnm_fc(fdata, Nrow, Ncol, 
                       Nmax, -m,
                       pt, Nmax - abs(m) + 1, 
                       s, Q, 
                       ff, Q, 
                       hkm, Nmax+1, 
                       y,Nmax+1,
                       kiss_cfg_fw,
                       kiss_cfg_bw);

		pt = sc + inds[2*m];
		bnm_fc(fdata, Nrow, Ncol,
                       Nmax, m,
                       pt, Nmax - abs(m) + 1,
                       s, Q,
                       ff, Q,
                       hkm, Nmax+1,
                       y, Nmax+1,
                       kiss_cfg_fw,
                       kiss_cfg_bw);
	}

        //Free all memory here
	free(s);
        free(kiss_cfg_fw);
        free(kiss_cfg_bw);
        free(hkm);
        free(y);
        free(ff);
        free(inds);
}

