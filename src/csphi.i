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

%apply (double* IN_ARRAY1, int DIM1) {(double* y, int len)};

%apply (SCOMPLEX* IN_ARRAY1, int DIM1) {(SCOMPLEX* s,int Q)};

%include "csphi.h"
 
