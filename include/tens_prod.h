#ifndef TENS_PROD_H
#define TENS_PROD_H

#include<complex.h>
#include<stdio.h>

#include"macro.h"

// see Luscher Weisz JHEP 0109 p. 010 (2001)   (hep-lat/0108014)

typedef struct TensProd {
   double complex comp[NCOLOR][NCOLOR][NCOLOR][NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
} TensProd;

void zero_TensProd(TensProd * A);
void one_TensProd(TensProd *A);

void equal_TensProd(TensProd * restrict A, TensProd const * restrict B); // A=B

void times_equal_real_TensProd(TensProd * restrict A, double r); // A*=r
void times_equal_complex_TensProd(TensProd * restrict A, double complex r); // A*=r

void plus_equal_TensProd(TensProd * restrict A, TensProd const * restrict B); // A+=B

void times_TensProd(TensProd * restrict A,
                    TensProd const * restrict B,
                    TensProd const * restrict C); // A=B*C
void times_equal_TensProd(TensProd * restrict A, TensProd const * restrict B); // A*=B

double retr_TensProd(TensProd const * restrict A);
double imtr_TensProd(TensProd const * restrict A);

void print_on_screen_TensProd(TensProd const * const A);
void print_on_file_TensProd(FILE *fp, TensProd const * const A);
void print_on_binary_file_noswap_TensProd(FILE *fp, TensProd const * const A);
void print_on_binary_file_swap_TensProd(FILE *fp, TensProd const * const A);
void print_on_binary_file_bigen_TensProd(FILE *fp, TensProd const * const A);
void read_from_file_TensProd(FILE *fp, TensProd *A);
void read_from_binary_file_noswap_TensProd(FILE *fp, TensProd *A);
void read_from_binary_file_swap_TensProd(FILE *fp, TensProd *A);
void read_from_binary_file_bigen_TensProd(FILE *fp, TensProd *A);

#endif // TENS_PROD_H

