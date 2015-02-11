/* csphi.i */
%module csphi
%{
#define SWIG_FILE_WITH_INIT
#include "csphi.h"
%}

%include "numpy.i"

%init %{
import_array();
%}

%numpy_typemaps(SCOMPLEX, NPY_CDOUBLE, int)

%apply (double* INPLACE_ARRAY1, int DIM1) {(double* y, int len)};

%apply (SCOMPLEX* INPLACE_ARRAY1, int DIM1) {(SCOMPLEX* s, int Q)};

%apply (SCOMPLEX* IN_ARRAY2, int DIM1, int DIM2) {(SCOMPLEX* fdata, int Nrow, int Ncol)};

%apply (SCOMPLEX* INPLACE_ARRAY1, int DIM1) {(SCOMPLEX* sc, int L)};

%include "csphi.h"
 
