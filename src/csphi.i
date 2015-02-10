/* csphi.i */
%module csphi
%{
#define SWIG_FILE_WITH_INIT
#include <complex.h>
#include "csphi.h"
%}

%include <complex.i>
%include "numpy.i"

%init %{
import_array();
%}

%numpy_typemaps(float complex, NPY_CFLOAT , int)
%numpy_typemaps(double complex, NPY_CDOUBLE, int)

%apply (double* IN_ARRAY1, int DIM1) {(double* y, int len)};

%apply (double complex* IN_ARRAY1, int DIM1) {(double complex* s,int Q)};

%include "csphi.h"
 
