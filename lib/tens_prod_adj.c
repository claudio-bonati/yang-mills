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
void zero_TensProd_adj(TensProdAdj * restrict A);


// initialize to one
void one_TensProd_adj(TensProdAdj * restrict A);


// A=B
void equal_TensProd_adj(TensProdAdj * restrict A, TensProdAdj const * const restrict B);


// A*=r real
void times_equal_real_TensProd_adj(TensProdAdj * restrict A, double r);

// A+=B
void plus_equal_TensProd_adj(TensProdAdj * restrict A, TensProdAdj const * const restrict B);


// A=B*C
void times_TensProd_adj(TensProdAdj * restrict A,
                        TensProdAdj const * const restrict B,
                        TensProdAdj const * const restrict C);


// A*=B
void times_equal_TensProd_adj(TensProdAdj * restrict A, TensProdAdj const * const restrict B);


double retr_TensProd_adj(TensProdAdj const * const restrict A);
double imtr_TensProd_adj(TensProdAdj const * const restrict A);

#endif
