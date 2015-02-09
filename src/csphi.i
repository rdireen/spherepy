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

%apply (double* IN_ARRAY1, int DIM1) {(double* y, int len)};
%include "csphi.h"
 
