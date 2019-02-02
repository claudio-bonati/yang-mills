#ifndef TENS_PROD_ADJ_C
#define TENS_PROD_ADJ_C

#include<complex.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/aligncheck.h"
#include"../include/endianness.h"
#include"../include/macro.h"
#include"../include/tens_prod_adj.h"


// initialize to zero
void zero_TensProdAdj(TensProdAdj * restrict A);


// initialize to one
void one_TensProdAdj(TensProdAdj * restrict A);


// A=B
void equal_TensProdAdj(TensProdAdj * restrict A, TensProdAdj const * const restrict B);


// A*=r real
void times_equal_real_TensProdAdj(TensProdAdj * restrict A, double r);

// A+=B
void plus_equal_TensProdAdj(TensProdAdj * restrict A, TensProdAdj const * const restrict B);


// A=B*C
void times_TensProdAdj(TensProdAdj * restrict A,
                        TensProdAdj const * const restrict B,
                        TensProdAdj const * const restrict C);


// A*=B
void times_equal_TensProdAdj(TensProdAdj * restrict A, TensProdAdj const * const restrict B);


double retr_TensProdAdj(TensProdAdj const * const restrict A);
double imtr_TensProdAdj(TensProdAdj const * const restrict A);

#endif
