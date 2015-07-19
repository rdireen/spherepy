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

/*****************************************************************************


                         csphi: Low Level Routines 

Randy Direen
2/11/2015

These routines are used in the process of calculating the scalar 
spherical harmonic coefficients from spherical pattern information.
For python versions, see pysphi.py.

 *****************************************************************************/

#ifndef __CSPHI__
#define __CSPHI__

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "kiss_fft.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SFLOAT double
#define SINT int

    typedef struct {
        SFLOAT r;
        SFLOAT i;
    } SCOMPLEX;

#define SUCCESS    0
#define OUTBOUNDS  1

#define PI 3.14159265358979323846

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
    SFLOAT ynnm(int n, int m);
    
    
    //Some helper functions to use pointers in python
    SFLOAT *new_SFLOAT(SFLOAT val);
    SINT *new_SINT(SINT val);
    SFLOAT get_SFLOAT(SFLOAT *d);
    SINT get_SINT(SINT *i);
    void delete_SFLOAT(SFLOAT *d);
    void delete_SINT(SINT *i);

    /** 
     * Calculate the starting value y(n|nm) for the recursion formula
     *that generates the Fourier coefficients of the spherical 
     *harmonic Ynm. This is the high dynamic range version used for calculating
     *modes larger than n = 1000.
     * @brief Calculate y(n|nm) 
     * 
     * @param n The n index in Ynm
     * @param m The m index in Ynm
     * @param val is the starting value y(n|nm), but needs to be normalized by ee
     * @param ee is the integer used to normalize y(n|nm) with 10^ee.
     * 
     * @see ynunm_hdr()
     */
    void ynnm_hdr(int n, int m, SFLOAT *val, SINT *ee);

    /** 
     * Calculates the Fourier coefficients y(nu|n,m). y must have a
     * length greater than n.  
     * @brief Calculate y(nu|nm) 
     * 
     * @param en Index n
     * @param em Index m
     * @param y  Array of Fourier coefficients for Ynm
     * @param len Length of array y
     * 
     * @see ynnm()
     */
    void ynunm(int en, int em, SFLOAT* y, int len);

    /** 
     * Calculates the Fourier coefficients y(nu|n,m). y must have a
     * length greater than n. This is the high dynamic range version that
     * works for n greater than 1000.  
     * @brief Calculate y(nu|nm) 
     * 
     * @param en Index n
     * @param em Index m
     * @param y  Array of Fourier coefficients for Ynm
     * @param len Length of array y
     * 
     * @see ynnm()
     */
    void ynunm_hdr(int en, int em, SINT* EE, int len_e, SFLOAT* y, int len);

    /** 
     * @brief Finds the next largest number to S, that can be factored into small
     * primes.
     * 
     * @param S The output will be a number greater than or equal to this
     * 
     * @return A number greater than S, that can be factored by 2,3,5
     */
    int FindQ(int S);

    /** 
     * @brief Calculate s data to be passed to hkm_fc
     * 
     * @param s  a vector of the s data with length Q
     * @param Q
     * @param Nrows number of rows in fdata
     * @param Nmax largest spherical coefficient of interest
     */
    void SData(SCOMPLEX* s, int Q, int Nrows, int NcoefMax);

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
     * @param hkm   Output 
     * @param len   Length of output, must be greater than n and less than or equal to Q (set to n+1)
     * @param ss    Fourier transformed version of the data provided by SData
     * @param Q     A number larger than Nrow+n, make factorable into small primes for speed
     * @param Q2    Same as Q (used this so the swig interface creator wouldn't get confused)
     * @param kiss_cfg_fw    Config info for a forwards FFT
     * @param kiss_cfg_bw    Config info for a backwards FFT
     */
    void hkm_fc(SCOMPLEX* gcoef, int Nrow, int Ncol,
            int n, int m,
            SCOMPLEX* hkm, int len,
            SCOMPLEX* ss, int Q,
            SCOMPLEX* ff, int Q2,
            kiss_fft_cfg kiss_cfg_fw,
            kiss_fft_cfg kiss_cfg_bw);

    /** 
     * @brief Calcutates a vector vec of spherical coefficients
     * 
     * @param fdata Fourier coefficient data
     * @param Nrow  Number of rows in the Fourier coefficient data
     * @param Ncol  Number of columns in the Fourier coefficient data
     * @param Nmax Maximum number of spherical coefficients desired
     * @param m  The m index into the spherical coefficients
     * @param vec The output vector of length Nmax-abs(m)+1
     * @param ss Fourier transform of the data returned from SData
     * @param Q The length of ss
     * @param ff Some work space that is as long as ss
     * @param Q2 Same as Q, I included it so swig could easily be used
     * @param hkm See hkm_fc
     * @param Lhkm Length of hkm
     * @param y See ynunm
     * @param Ly Length of y
     * @param kiss_cfg_fw Config info for a forwards FFT
     * @param kiss_cfg_bw Config info for a backwards FFT
     */
    void bnm_fc(SCOMPLEX * fdata, int Nrow, int Ncol,
            int Nmax, int m,
            SCOMPLEX* vec, int L,
            SCOMPLEX* ss, int Q,
            SCOMPLEX* ff, int Q2,
            SCOMPLEX* hkm, int Lhkm,
            SFLOAT* y, int Ly,
            kiss_fft_cfg kiss_cfg_fw,
            kiss_fft_cfg kiss_cfg_bw);

    /** 
     * This the high dynamic range version of bnm_fc that works for n > 1000
     * @brief Calcutates a vector vec of spherical coefficients
     * 
     * @param fdata Fourier coefficient data
     * @param Nrow  Number of rows in the Fourier coefficient data
     * @param Ncol  Number of columns in the Fourier coefficient data
     * @param Nmax Maximum number of spherical coefficients desired
     * @param m  The m index into the spherical coefficients
     * @param vec The output vector of length Nmax-abs(m)+1
     * @param ss Fourier transform of the data returned from SData
     * @param Q The length of ss
     * @param ff Some work space that is as long as ss
     * @param Q2 Same as Q, I included it so swig could easily be used
     * @param hkm See hkm_fc
     * @param Lhkm Length of hkm
     * @param y See ynunm
     * @param Ly Length of y
     * @param kiss_cfg_fw Config info for a forwards FFT
     * @param kiss_cfg_bw Config info for a backwards FFT
     */
    void bnm_fc_hdr(SCOMPLEX * fdata, int Nrow, int Ncol,
            int Nmax, int m,
            SCOMPLEX* vec, int L,
            SCOMPLEX* ss, int Q,
            SCOMPLEX* ff, int Q2,
            SCOMPLEX* hkm, int Lhkm,
            SINT* EE, int Le,
            SFLOAT* y, int Ly,
            kiss_fft_cfg kiss_cfg_fw,
            kiss_fft_cfg kiss_cfg_bw);

    /** 
     * @brief Fourier coefficients to spherical harmonic coefficients
     * 
     * @param fdata Fourier coefficients
     * @param Nrow Number of rows in fdata
     * @param Ncol Number of columns in fdata
     * @param sc Array of coefficients
     * @param L Length of sc
     * @param Nmax Maximum number of n modes desired
     * @param Mmax Maximum number of abs(m) modes desired
     */
    void fc_to_sc(SCOMPLEX* fdata, int Nrow, int Ncol,
            SCOMPLEX* sc, int L,
            int Nmax, int Mmax);

    /** 
     * This is the high dynamic range version that works for n > 1000
     * @brief Fourier coefficients to spherical harmonic coefficients
     * 
     * @param fdata Fourier coefficients
     * @param Nrow Number of rows in fdata
     * @param Ncol Number of columns in fdata
     * @param sc Array of coefficients
     * @param L Length of sc
     * @param Nmax Maximum number of n modes desired
     * @param Mmax Maximum number of abs(m) modes desired
     */
    void fc_to_sc_hdr(SCOMPLEX* fdata, int Nrow, int Ncol,
            SCOMPLEX* sc, int L,
            int Nmax, int Mmax);

    /** 
     * @brief Fills one of the columns of the Fourier coefficient array using the spherical harmonic coefficients
     * 
     * @param vec Array of spherical harmonic coefficients
     * @param m The m index
     * @param Nmax Largest spherical harmonic coefficient
     * @param fdata Fourier coefficients
     * @param Nrow Number of rows in the Fourier coefficients array
     * @param Ncol Number of columns in the Fourier coefficients array
     * @param M m index into the Fourier coefficients array
     * @param y output of ynunm
     * @param len length of y
     */
    void fcvec_m_sc(SCOMPLEX* vec,
            int m, int Nmax,
            SCOMPLEX* fdata, int Nrow, int Ncol,
            int M,
            SFLOAT * y, int len);

    /** 
     * This is the high dynamic range version that works for n > 1000
     * @brief Fills one of the columns of the Fourier coefficient array using the spherical harmonic coefficients
     * 
     * @param vec Array of spherical harmonic coefficients
     * @param m The m index
     * @param Nmax Largest spherical harmonic coefficient
     * @param fdata Fourier coefficients
     * @param Nrow Number of rows in the Fourier coefficients array
     * @param Ncol Number of columns in the Fourier coefficients array
     * @param M m index into the Fourier coefficients array
     * @param y output of ynunm
     * @param len length of y
     */
    void fcvec_m_sc_hdr(SCOMPLEX* vec,
            int m, int Nmax,
            SCOMPLEX* fdata, int Nrow, int Ncol,
            int M,
            SINT * EE, int len_e,
            SFLOAT * y, int len);

    /** 
     * @brief Spherical harmonic coefficients to Fourier coefficients
     * 
     * @param fdata Fourier coefficients
     * @param Nrow Number of rows in fdata
     * @param Ncol Number of columns in fdata
     * @param sc Array of coefficients
     * @param L Length of sc
     * @param Nmax Maximum number of n modes in sc
     * @param Mmax Maximum number of abs(m) modes in sc
     */
    void sc_to_fc(SCOMPLEX* fdata, int Nrow, int Ncol,
            SCOMPLEX* sc, int L,
            int Nmax, int Mmax);

    /** 
     * This is the high dynamic range version that works for n > 1000
     * @brief Spherical harmonic coefficients to Fourier coefficients
     * 
     * @param fdata Fourier coefficients
     * @param Nrow Number of rows in fdata
     * @param Ncol Number of columns in fdata
     * @param sc Array of coefficients
     * @param L Length of sc
     * @param Nmax Maximum number of n modes in sc
     * @param Mmax Maximum number of abs(m) modes in sc
     */
    void sc_to_fc_hdr(SCOMPLEX* fdata, int Nrow, int Ncol,
            SCOMPLEX* sc, int L,
            int Nmax, int Mmax);

#ifdef __cplusplus
}
#endif

#endif //__CSPHI__
