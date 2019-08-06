#ifndef FLAVOUR_MATRIX_C
#define FLAVOUR_MATRIX_C

#include<complex.h>
#include<math.h>

#include"../include/macro.h"
#include"../include/flavour_matrix.h"

// A=0
void zero_FMatrix(FMatrix * restrict A);


// A=B
void equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B);


// A+=B
void plus_equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B);


// A-=B
void minus_equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B);


// A-=B^{dag}
void minus_equal_dag_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B);


// A*=r
void times_equal_real_FMatrix(FMatrix * restrict A, double r);


// A*=r
void times_equal_complex_FMatrix(FMatrix * restrict A, double complex r);


// A*=B
void times_equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B);


// real part of the trace /NHIGGS
double retr_FMatrix(FMatrix const * const restrict A);


// imaginary part of the trace /NHIGGS
double imtr_FMatrix(FMatrix const * const restrict A);


// l2 norm of the matrix
double norm_FMatrix(FMatrix const * const restrict A);



#endif

