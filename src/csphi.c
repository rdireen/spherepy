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

#include "csphi.h"

/**
 * This is a MACRO that helps in the ynunm routine
 * @brief MACRO for ynunm routine
 * @see ynunm()
 */
#define B(n,v) (((n)-(v))*((n)+(v)+1.0))

/*----------------------------------------------------------------------------*/
double ynnm(int n, int m) {
    int pm, k;
    double a, ynnm;

    a = 1.0 / sqrt(4.0 * PI);
    pm = abs(m);

    if (n < pm) {
        ynnm = 0;
    } else if (n == 0) {
        ynnm = a;
    } else {
        ynnm = a;
        for (k = 1; k <= n; k++) {
            ynnm = sqrt((2.0 * k + 1) / (8.0 * k)) * ynnm;
        }
        if (n != pm) {
            for (k = n - 1; k >= pm; k--) {
                ynnm = sqrt(1.0 * (n + k + 1.0) / (n - k)) * ynnm;
            }
        }
    }
    return ynnm;
}

/*----------------------------------------------------------------------------*/
SFLOAT *new_SFLOAT(SFLOAT val) {
    SFLOAT *d = (SFLOAT *) malloc(sizeof (val));
    *d = val;
    return d;
}

SINT *new_SINT(SINT val) {
    SINT *i = (SINT *) malloc(sizeof (val));
    *i = val;
    return i;
}

SFLOAT get_SFLOAT(SFLOAT *d) {
    return *d;
}

SINT get_SINT(SINT *i) {
    return *i;
}

void delete_SFLOAT(SFLOAT *d) {
    free(d);
}

void delete_SINT(SINT *i) {
    free(i);
}

void ynnm_hdr(int n, int m, SFLOAT *val, SINT *ee) {
    int pm, k;
    double a, ynnm;
    int cst;
    double low_bound, high_bound;

    a = 1.0 / sqrt(4.0 * PI);
    pm = abs(m);

    cst = 100;
    (*ee) = 0;
    low_bound = 1e-100;
    high_bound = 1e100;

    if (n < pm) {
        ynnm = 0;
    } else if (n == 0) {
        ynnm = a;
    } else {
        ynnm = a;
        for (k = 1; k <= n; k++) {
            ynnm = sqrt((2.0 * k + 1) / (8.0 * k)) * ynnm;
            if (ynnm < low_bound) {
                ynnm *= pow(10, cst);
                (*ee) += cst;
            }
        }
        if (n != pm) {
            for (k = n - 1; k >= pm; k--) {
                ynnm = sqrt(1.0 * (n + k + 1.0) / (n - k)) * ynnm;
                if (ynnm > high_bound) {
                    ynnm *= pow(10, -cst);
                    (*ee) -= cst;
                }
            }
        }
    }
    (*val) = ynnm;
}

/*----------------------------------------------------------------------------*/
void ynunm(int en, int em, SFLOAT* y, int len) {
    int k;

    for (k = 0; k < len; k++)
        *(y + k) = 0;

    if (abs(em) <= en) {
        *(y + en) = ynnm(en, em);
        k = en - 2;
        if (k >= 0) {
            *(y + k) = (B(en, k + 1.0) + B(en, k + 2.0) - 4.0 * em * em)*(*(y + k + 2)) / B(en, k);
            for (k = k - 2; k >= 0; k -= 2) {
                (*(y + k)) = ((B(en, k + 1.0) + B(en, k + 2.0) - 4.0 * em * em)*(*(y + k + 2)) - B(en, k + 3.0)*(*(y + k + 4))) / B(en, k);
            }
        }
    }
}

/*----------------------------------------------------------------------------*/
void ynunm_hdr(int en, int em, SINT* EE, int len_e, SFLOAT* y, int len) {

    SFLOAT ynnm_norm = 0;
    SINT e1 = 0;
    SFLOAT *val_p = &ynnm_norm;
    SINT *ee_p = &e1;

    int cst, ee;
    double high_bound;

    double exp10 = 1.0;
    int mod_odd = 0;

    int k;

    cst = 100;
    ee = 0;
    high_bound = 1e100;

    ynnm_hdr(en, em, val_p, ee_p);


    for (k = 0; k < len; k++)
        *(y + k) = 0;

    if (abs(em) <= en) {
        *(y + en) = 1.0; // Start from 1.0 instead of ynnm 
        k = en - 2;
        if (k >= 0) {
            *(y + k) = (B(en, k + 1.0) + B(en, k + 2.0) - 4.0 * em * em)*(*(y + k + 2)) / B(en, k);
            for (k = k - 2; k >= 0; k -= 2) {
                (*(y + k)) = ((B(en, k + 1.0) + B(en, k + 2.0) - 4.0 * em * em)*(*(y + k + 2)) - B(en, k + 3.0)*(*(y + k + 4))) / B(en, k);
                *(EE + k) = ee;
                if (*(y + k) > high_bound) {
                    *(y + k) *= pow(10, -cst);
                    *(y + k + 2) *= pow(10, -cst);
                    ee -= cst;
                    *(EE + k) = ee;
                    *(EE + k + 2) = ee;
                }
            }
        }
    }

    mod_odd = en % 2;
    for (k = mod_odd; k < len; k += 2) {
        exp10 = pow(10, *(EE + k) + e1);
        *(y + k) *= ynnm_norm / exp10;
    }
}

/*----------------------------------------------------------------------------*/
int FindQ(int S) {
    int A;
    int Q = S;

    A = Q;
    while (A != 1) {
        if (A % 2 == 0)
            A = A / 2;
        else if (A % 3 == 0)
            A = A / 3;
        else if (A % 5 == 0)
            A = A / 5;
        else {
            A = Q + 1;
            Q = A;
        }
    }

    return Q;
}

/*----------------------------------------------------------------------------*/
void SData(SCOMPLEX* s, int Q, int Nrows, int Nmax) {
    int mm, k;
    int nn = Nmax + 1;

    memset(s, 0, Q * sizeof (SCOMPLEX));

    //we assume Nrows is even
    mm = Nrows / 2;

    for (k = mm; k <= mm + nn - 1; k++)
        if ((k % 2) == 1)
            s[k - mm].i = -1 / ((SFLOAT) k);
    for (k = -mm; k < mm; k++)
        if ((abs(k) % 2) == 1)
            s[Q + k - mm].i = -1 / ((SFLOAT) k);
}

/*----------------------------------------------------------------------------*/
void hkm_fc(SCOMPLEX* gcoef, int Nrow, int Ncol,
        int n, int m,
        SCOMPLEX* hkm, int len,
        SCOMPLEX* ss, int Q,
        SCOMPLEX* ff, int Q2,
        kiss_fft_cfg kiss_cfg_fw,
        kiss_fft_cfg kiss_cfg_bw) {
    int k, ind, mm;
    SFLOAT ffr;

    memset(ff, 0, Q * sizeof (SCOMPLEX));
    memset(hkm, 0, len * sizeof (SCOMPLEX));

    // properly get the index into gcoef that m corresponds to
    if (m >= 0)
        ind = m;
    else
        ind = Ncol + m;

    // we assume the number of rows is even
    mm = Nrow / 2;

    for (k = 0; k < mm; k++) {
        ff[k] = gcoef[Ncol * (mm + k) + ind];
        ff[mm + k] = gcoef[Ncol * k + ind];
    }

    kiss_fft(kiss_cfg_fw, (kiss_fft_cpx*) ff, (kiss_fft_cpx*) ff);

    for (k = 0; k < Q; k++) {
        ffr = ff[k].r;
        ff[k].r = ss[k].r * ff[k].r - ss[k].i * ff[k].i;
        ff[k].i = ss[k].r * ff[k].i + ss[k].i*ffr;
    }

    kiss_fft(kiss_cfg_bw, (kiss_fft_cpx*) ff, (kiss_fft_cpx*) ff);

    for (k = 0; k < len; k++) {
        hkm[k].r = 4.0 * PI * ff[k].r / (SFLOAT) Q;
        hkm[k].i = 4.0 * PI * ff[k].i / (SFLOAT) Q;
    }

}

/*----------------------------------------------------------------------------*/
void bnm_fc(SCOMPLEX * fdata, int Nrow, int Ncol,
        int Nmax, int m,
        SCOMPLEX* vec, int L,
        SCOMPLEX* ss, int Q,
        SCOMPLEX* ff, int Q2,
        SCOMPLEX* hkm, int Lhkm,
        SFLOAT* y, int Ly,
        kiss_fft_cfg kiss_cfg_fw,
        kiss_fft_cfg kiss_cfg_bw) {
    int n, k;
    int absm = abs(m);
    SFLOAT vecr;

    hkm_fc(fdata, Nrow, Ncol,
            Nmax, m,
            hkm, Lhkm,
            ss, Q,
            ff, Q2,
            kiss_cfg_fw,
            kiss_cfg_bw);

    for (n = absm; n <= Nmax; n++) {
        ynunm(n, m, y, Ly);

        vec[n - absm].r = hkm[0].r * y[0];
        vec[n - absm].i = hkm[0].i * y[0];

        for (k = 1; k < n + 1; k++) {
            //TODO: need to sort before I sum
            vec[n - absm].r += 2.0 * hkm[k].r * y[k];
            vec[n - absm].i += 2.0 * hkm[k].i * y[k];
        }

        //In the following: vec[n-absm] = i**(-m) * vec[n-absm] 
        if (absm % 2 == 1) {
            vecr = vec[n - absm].r;
            vec[n - absm].r = pow(-1.0, (m - 1) / 2) * vec[n - absm].i;
            vec[n - absm].i = -pow(-1.0, (m - 1) / 2) * vecr;
        } else {
            vec[n - absm].r = pow(-1.0, m / 2) * vec[n - absm].r;
            vec[n - absm].i = pow(-1.0, m / 2) * vec[n - absm].i;
        }
    }
}

/*----------------------------------------------------------------------------*/
void bnm_fc_hdr(SCOMPLEX * fdata, int Nrow, int Ncol,
        int Nmax, int m,
        SCOMPLEX* vec, int L,
        SCOMPLEX* ss, int Q,
        SCOMPLEX* ff, int Q2,
        SCOMPLEX* hkm, int Lhkm,
        SINT* EE, int Le,
        SFLOAT* y, int Ly,
        kiss_fft_cfg kiss_cfg_fw,
        kiss_fft_cfg kiss_cfg_bw) {
    int n, k;
    int absm = abs(m);
    SFLOAT vecr;

    hkm_fc(fdata, Nrow, Ncol,
            Nmax, m,
            hkm, Lhkm,
            ss, Q,
            ff, Q2,
            kiss_cfg_fw,
            kiss_cfg_bw);

    for (n = absm; n <= Nmax; n++) {
        memset(EE, 0, Le * sizeof (SINT));
        ynunm_hdr(n, m, EE, Le, y, Ly);

        vec[n - absm].r = hkm[0].r * y[0];
        vec[n - absm].i = hkm[0].i * y[0];

        for (k = 1; k < n + 1; k++) {
            //TODO: need to sort before I sum
            vec[n - absm].r += 2.0 * hkm[k].r * y[k];
            vec[n - absm].i += 2.0 * hkm[k].i * y[k];
        }

        //In the following: vec[n-absm] = i**(-m) * vec[n-absm] 
        if (absm % 2 == 1) {
            vecr = vec[n - absm].r;
            vec[n - absm].r = pow(-1.0, (m - 1) / 2) * vec[n - absm].i;
            vec[n - absm].i = -pow(-1.0, (m - 1) / 2) * vecr;
        } else {
            vec[n - absm].r = pow(-1.0, m / 2) * vec[n - absm].r;
            vec[n - absm].i = pow(-1.0, m / 2) * vec[n - absm].i;
        }
    }
}

/*----------------------------------------------------------------------------*/
void fc_to_sc(SCOMPLEX* fdata, int Nrow, int Ncol,
        SCOMPLEX* sc, int L,
        int Nmax, int Mmax) {
    int m, Q;
    SCOMPLEX* pt;
    int* inds;
    int N = Nmax + 1;
    int QQ = Nmax;
    SCOMPLEX* s;
    SCOMPLEX* hkm;
    SCOMPLEX* ff;
    SFLOAT* y;
    kiss_fft_cfg kiss_cfg_fw;
    kiss_fft_cfg kiss_cfg_bw;

    Q = FindQ(Nrow + Nmax);

    //Allocate all memory here
    s = (SCOMPLEX*) malloc(Q * sizeof (SCOMPLEX));
    SData(s, Q, Nrow, Nmax);

    kiss_cfg_fw = kiss_fft_alloc(Q, 0, 0, 0);
    kiss_cfg_bw = kiss_fft_alloc(Q, 1, 0, 0);

    kiss_fft(kiss_cfg_fw, (kiss_fft_cpx*) s, (kiss_fft_cpx*) s);

    hkm = (SCOMPLEX*) malloc((Nmax + 1) * sizeof (SCOMPLEX));
    y = (SFLOAT*) malloc((Nmax + 1) * sizeof (SFLOAT));
    ff = (SCOMPLEX*) malloc(Q * sizeof (SCOMPLEX));

    inds = (int*) malloc((2 * Mmax + 1) * sizeof (int));
    memset(inds, 0, (2 * Mmax + 1) * sizeof (int));


    //setup indices that index into sc
    inds[0] = 0;
    for (m = 1; m <= Mmax; m++) {
        inds[2 * m - 1] = N;
        N = N + QQ;
        inds[2 * m] = N;
        N = N + QQ;
        QQ--;
    }


    //Calculate each of the bnm coefficients
    bnm_fc(fdata, Nrow, Ncol,
            Nmax, 0,
            sc, Nmax + 1,
            s, Q,
            ff, Q,
            hkm, Nmax + 1,
            y, Nmax + 1,
            kiss_cfg_fw,
            kiss_cfg_bw);

    for (m = 1; m <= Mmax; m++) {
        pt = sc + inds[2 * m - 1];
        bnm_fc(fdata, Nrow, Ncol,
                Nmax, -m,
                pt, Nmax - abs(m) + 1,
                s, Q,
                ff, Q,
                hkm, Nmax + 1,
                y, Nmax + 1,
                kiss_cfg_fw,
                kiss_cfg_bw);

        pt = sc + inds[2 * m];
        bnm_fc(fdata, Nrow, Ncol,
                Nmax, m,
                pt, Nmax - abs(m) + 1,
                s, Q,
                ff, Q,
                hkm, Nmax + 1,
                y, Nmax + 1,
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

/*----------------------------------------------------------------------------*/
void fc_to_sc_hdr(SCOMPLEX* fdata, int Nrow, int Ncol,
        SCOMPLEX* sc, int L,
        int Nmax, int Mmax) {
    int m, Q;
    SCOMPLEX* pt;
    int* inds;
    int N = Nmax + 1;
    int QQ = Nmax;
    SCOMPLEX* s;
    SCOMPLEX* hkm;
    SCOMPLEX* ff;
    SINT* EE;
    SFLOAT* y;
    kiss_fft_cfg kiss_cfg_fw;
    kiss_fft_cfg kiss_cfg_bw;

    Q = FindQ(Nrow + Nmax);

    //Allocate all memory here
    s = (SCOMPLEX*) malloc(Q * sizeof (SCOMPLEX));
    SData(s, Q, Nrow, Nmax);

    kiss_cfg_fw = kiss_fft_alloc(Q, 0, 0, 0);
    kiss_cfg_bw = kiss_fft_alloc(Q, 1, 0, 0);

    kiss_fft(kiss_cfg_fw, (kiss_fft_cpx*) s, (kiss_fft_cpx*) s);

    hkm = (SCOMPLEX*) malloc((Nmax + 1) * sizeof (SCOMPLEX));
    EE = (SINT*) malloc((Nmax + 1) * sizeof (SINT));
    y = (SFLOAT*) malloc((Nmax + 1) * sizeof (SFLOAT));
    ff = (SCOMPLEX*) malloc(Q * sizeof (SCOMPLEX));

    inds = (int*) malloc((2 * Mmax + 1) * sizeof (int));
    memset(inds, 0, (2 * Mmax + 1) * sizeof (int));


    //setup indices that index into sc
    inds[0] = 0;
    for (m = 1; m <= Mmax; m++) {
        inds[2 * m - 1] = N;
        N = N + QQ;
        inds[2 * m] = N;
        N = N + QQ;
        QQ--;
    }

    //Calculate each of the bnm coefficients
    bnm_fc_hdr(fdata, Nrow, Ncol,
            Nmax, 0,
            sc, Nmax + 1,
            s, Q,
            ff, Q,
            hkm, Nmax + 1,
            EE, Nmax + 1,
            y, Nmax + 1,
            kiss_cfg_fw,
            kiss_cfg_bw);

    for (m = 1; m <= Mmax; m++) {
        pt = sc + inds[2 * m - 1];
        bnm_fc_hdr(fdata, Nrow, Ncol,
                Nmax, -m,
                pt, Nmax - abs(m) + 1,
                s, Q,
                ff, Q,
                hkm, Nmax + 1,
                EE, Nmax + 1,
                y, Nmax + 1,
                kiss_cfg_fw,
                kiss_cfg_bw);

        pt = sc + inds[2 * m];
        bnm_fc_hdr(fdata, Nrow, Ncol,
                Nmax, m,
                pt, Nmax - abs(m) + 1,
                s, Q,
                ff, Q,
                hkm, Nmax + 1,
                EE, Nmax + 1,
                y, Nmax + 1,
                kiss_cfg_fw,
                kiss_cfg_bw);
    }

    //Free all memory here
    free(s);
    free(kiss_cfg_fw);
    free(kiss_cfg_bw);
    free(hkm);
    free(EE);
    free(y);
    free(ff);
    free(inds);
}

/*----------------------------------------------------------------------------*/
void fcvec_m_sc(SCOMPLEX * vec,
        int m, int Nmax,
        SCOMPLEX * fdata, int Nrow, int Ncol,
        int M,
        SFLOAT* y, int len) {

    int n, k, mm, H;
    int K = Nmax + 1;
    int absm = abs(m);
    SCOMPLEX * pfdata;
    SCOMPLEX tmp;
    SFLOAT tmpr;

    for (n = absm; n <= Nmax; n++) {
        ynunm(n, m, y, len);

        //TODO: this loop can be made faster
        for (k = 0; k < K; k++) {
            pfdata = fdata + k * Ncol + M;
            tmp = *pfdata;
            tmp.r += vec[n - absm].r * y[k];
            tmp.i += vec[n - absm].i * y[k];
            *pfdata = tmp;
        }
    }

    for (k = 0; k < K; k++) {
        pfdata = fdata + k * Ncol + M;
        tmp = *pfdata;

        //In the following: tmp = i**(m) * tmp 
        if (abs(m) % 2 == 1) {
            tmpr = tmp.r;
            tmp.r = -pow(-1.0, (m - 1) / 2) * tmp.i;
            tmp.i = pow(-1.0, (m - 1) / 2) * tmpr;
        } else {
            tmp.r = pow(-1.0, m / 2) * tmp.r;
            tmp.i = pow(-1.0, m / 2) * tmp.i;
        }

        *pfdata = tmp;
    }


    //apply symmetry
    mm = pow(-1, m);

    if ((Nrow % 2) == 0)
        H = Nrow / 2 - 1;
    else
        H = (Nrow - 1) / 2;

    for (k = 0; k < H; k++) {
        fdata[Ncol * (Nrow - 1 - k) + M].r = mm * fdata[Ncol * (k + 1) + M].r;
        fdata[Ncol * (Nrow - 1 - k) + M].i = mm * fdata[Ncol * (k + 1) + M].i;
    }
}

/*----------------------------------------------------------------------------*/
void fcvec_m_sc_hdr(SCOMPLEX * vec,
        int m, int Nmax,
        SCOMPLEX * fdata, int Nrow, int Ncol,
        int M,
        SINT* EE, int len_e,
        SFLOAT* y, int len) {
    int n, k, mm, H;
    int K = Nmax + 1;
    int absm = abs(m);
    SCOMPLEX * pfdata;
    SCOMPLEX tmp;
    SFLOAT tmpr;

    for (n = absm; n <= Nmax; n++) {
        memset(EE, 0, len_e * sizeof (SINT));
        ynunm_hdr(n, m, EE, len_e, y, len);

        //TODO: this loop can be made faster
        for (k = 0; k < K; k++) {
            pfdata = fdata + k * Ncol + M;
            tmp = *pfdata;
            tmp.r += vec[n - absm].r * y[k];
            tmp.i += vec[n - absm].i * y[k];
            *pfdata = tmp;
        }
    }

    for (k = 0; k < K; k++) {
        pfdata = fdata + k * Ncol + M;
        tmp = *pfdata;

        //In the following: tmp = i**(m) * tmp 
        if (abs(m) % 2 == 1) {
            tmpr = tmp.r;
            tmp.r = -pow(-1.0, (m - 1) / 2) * tmp.i;
            tmp.i = pow(-1.0, (m - 1) / 2) * tmpr;
        } else {
            tmp.r = pow(-1.0, m / 2) * tmp.r;
            tmp.i = pow(-1.0, m / 2) * tmp.i;
        }

        *pfdata = tmp;
    }


    //apply symmetry
    mm = pow(-1, m);

    if ((Nrow % 2) == 0)
        H = Nrow / 2 - 1;
    else
        H = (Nrow - 1) / 2;

    for (k = 0; k < H; k++) {
        fdata[Ncol * (Nrow - 1 - k) + M].r = mm * fdata[Ncol * (k + 1) + M].r;
        fdata[Ncol * (Nrow - 1 - k) + M].i = mm * fdata[Ncol * (k + 1) + M].i;
    }
}

/*----------------------------------------------------------------------------*/
void sc_to_fc(SCOMPLEX* fdata, int Nrow, int Ncol,
        SCOMPLEX* sc, int L,
        int Nmax, int Mmax) {
    int k, m;
    int* inds;
    SFLOAT* y;
    int N = Nmax + 1;
    int QQ = Nmax;


    y = (SFLOAT*) malloc((Nmax + 1) * sizeof (SFLOAT));
    inds = (int*) malloc((2 * Mmax + 1) * sizeof (int));
    memset(inds, 0, (2 * Mmax + 1) * sizeof (int));

    //setup indices that index into sc
    inds[0] = 0;
    for (m = 1; m <= Mmax; m++) {
        inds[2 * m - 1] = N;
        N = N + QQ;
        inds[2 * m] = N;
        N = N + QQ;
        QQ--;
    }


    for (k = 0; k <= floor(Ncol / 2); k++) {
        if (k < Mmax) {
            fcvec_m_sc(sc + inds[2 * k],
                    k, Nmax,
                    fdata, Nrow, Ncol,
                    k,
                    y, Nmax + 1);
            fcvec_m_sc(sc + inds[2 * (k + 1) - 1],
                    -(k + 1), Nmax,
                    fdata, Nrow, Ncol,
                    Ncol - 1 - k,
                    y, Nmax + 1);
        } else if (k == Mmax) {
            fcvec_m_sc(sc + inds[2 * k],
                    k, Nmax,
                    fdata, Nrow, Ncol,
                    k,
                    y, Nmax + 1);
        }
    }


    free(y);
    free(inds);
}

/*----------------------------------------------------------------------------*/
void sc_to_fc_hdr(SCOMPLEX* fdata, int Nrow, int Ncol,
        SCOMPLEX* sc, int L,
        int Nmax, int Mmax) {
    int k, m;
    int* inds;
    SINT* EE;
    SFLOAT* y;
    int N = Nmax + 1;
    int QQ = Nmax;

    EE = (SINT*) malloc((Nmax + 1) * sizeof (SINT));
    y = (SFLOAT*) malloc((Nmax + 1) * sizeof (SFLOAT));
    inds = (int*) malloc((2 * Mmax + 1) * sizeof (int));
    memset(inds, 0, (2 * Mmax + 1) * sizeof (int));

    //setup indices that index into sc
    inds[0] = 0;
    for (m = 1; m <= Mmax; m++) {
        inds[2 * m - 1] = N;
        N = N + QQ;
        inds[2 * m] = N;
        N = N + QQ;
        QQ--;
    }


    for (k = 0; k <= floor(Ncol / 2); k++) {
        if (k < Mmax) {
            fcvec_m_sc_hdr(sc + inds[2 * k],
                    k, Nmax,
                    fdata, Nrow, Ncol,
                    k,
                    EE, Nmax + 1,
                    y, Nmax + 1);
            fcvec_m_sc_hdr(sc + inds[2 * (k + 1) - 1],
                    -(k + 1), Nmax,
                    fdata, Nrow, Ncol,
                    Ncol - 1 - k,
                    EE, Nmax + 1,
                    y, Nmax + 1);
        } else if (k == Mmax) {
            fcvec_m_sc_hdr(sc + inds[2 * k],
                    k, Nmax,
                    fdata, Nrow, Ncol,
                    k,
                    EE, Nmax + 1,
                    y, Nmax + 1);
        }
    }

    free(EE);
    free(y);
    free(inds);
}

/*----------------------------------------------------------------------------*/
int mindx(int m, int Nmax, int Mmax) {

    int ind = 0;
    int NN = Nmax + 1;
    int ii = 0;
    int i;

    if (m != 0) {
        ind = NN;
        ii = 1;
        for (i = 1; i < abs(m); i++) {
            ind = ind + 2 * (NN - i);
            ii = i + 1;
        }

        if (m > 0)
            ind = ind + NN - ii;
    }
    return ind;
}

/*----------------------------------------------------------------------------*/
void mode_nmajor_to_mmajor(int Nmax, int Mmax,
        SCOMPLEX* vec_nmajor, int len_nmajor,
        SCOMPLEX* vec_mmajor, int len_mmajor) {

    int n, m;
    int idx;
    int lim;
    int ii = 0;

    for (n = 0; n <= Nmax; n++) {
        if (n > Mmax) {
            lim = Mmax;
        } else {
            lim = n;
        }

        for (m = -lim; m <= lim; m++) {
            idx = mindx(m, Nmax, Mmax);
            *(vec_mmajor + ii + lim + m) = *(vec_nmajor + idx + n - abs(m));
        }

        ii += 2 * lim + 1;
    }
}

void mag_square_vec(SCOMPLEX* sc, int L, SFLOAT* out, int Lout) {

    int n = 0;
    SCOMPLEX* p;

    for (n = 0; n < L; n++) {
        p = sc + n;
        *(out + n) = (p->r * p->r) + (p->i * p->i);
    }
}

void abs_vec(SCOMPLEX* sc, int L, SFLOAT* out, int Lout) {

    int n = 0;
    SFLOAT* p;

    mag_square_vec(sc, L, out, Lout);

    for (n = 0; n < L; n++) {
        p = out + n;
        *(p) = sqrt(*p);
    }
}

int _cmp_func(const void * a, const void * b) {
    double aa = *(double*)a;
    double bb = *(double*)b;
    if (aa < bb) return -1;
    if (aa > bb) return 1;
    return 0;
}

void power_n(int Nmax, int Mmax, SCOMPLEX* sc, int L, SFLOAT* out, int Lout) {

    int n, m;
    int ii = 0;
    int lim;
    SFLOAT* work;

    work = (SFLOAT*) malloc(L * sizeof (SFLOAT));
    mag_square_vec(sc, L, work, L);

    for (n = 0; n <= Nmax; n++) {
        if (n > Mmax) {
            lim = Mmax;
        } else {
            lim = n;
        }

        if (lim == 0) {
            *(out + n) = sqrt(*(work + n));
        } else {
            qsort(work + ii, 2 * lim + 1, sizeof (SFLOAT), _cmp_func);
            for (m = -lim; m <= lim; m++) {
                *(out + n) += *(work + ii + lim + m);         
            }
        }
        ii += 2 * lim + 1;
    }

    free(work);
}

