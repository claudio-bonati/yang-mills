#ifndef U1_H
#define U1_H

#include<complex.h>
#include<stdio.h>

#include"macro.h"

typedef struct U1 {
     double complex comp __attribute__((aligned(DOUBLE_ALIGN)));
} U1;

void init_U1(U1 * restrict A, double complex vec);
void one_U1(U1 * restrict A);         // A=1
void zero_U1(U1 * restrict A);        // A=0

void equal_U1(U1 * restrict A, U1 const * restrict B);      // A=B
void equal_dag_U1(U1 * restrict A, U1 const * restrict B);  // A=B^{dag}

void plus_equal_U1(U1 * restrict A, U1 const * restrict B);     // A+=B
void plus_equal_dag_U1(U1 * restrict A, U1 const * restrict B); // A+=B^{dag}

void minus_equal_U1(U1 * restrict A, U1 const * restrict B);     // A-=B
void minus_equal_times_real_U1(U1 * restrict A, U1 const * restrict B, double r);     // A-=(r*B)
void minus_equal_dag_U1(U1 * restrict A, U1 const * restrict B); // A-=B^{dag}

void lin_comb_U1(U1 * restrict A,
                 double b, U1 const * restrict B,
                 double c, U1 const * restrict C);       // A=b*B+c*C
void lin_comb_dag1_U1(U1 * restrict A,
                      double b, U1 const * restrict B,
                      double c, U1 const * restrict C);  // A=b*B^{dag}+c*C
void lin_comb_dag2_U1(U1 * restrict A,
                      double b, U1 const * restrict B,
                      double c, U1 const * restrict C);  // A=b*B+c*C^{dag}
void lin_comb_dag12_U1(U1 *A,
                       double b, U1 const * restrict B,
                       double c, U1 const * restrict C); // A=b*B^{dag}+c*C^{dag}

void times_equal_real_U1(U1 * restrict A, double r); // A*=r

void times_equal_U1(U1 * restrict A, U1 const * restrict B);     // A*=B
void times_equal_dag_U1(U1 * restrict A, U1 const * restrict B); // A*=B^{dag}

void times_U1(U1 * restrict A,
              U1 const * restrict B,
              U1 const * restrict C);       // A=B*C
void times_dag1_U1(U1 * restrict A,
                   U1 const * restrict B,
                   U1 const * restrict C);  // A=B^{dag}*C
void times_dag2_U1(U1 * restrict A,
                   U1 const * restrict B,
                   U1 const * restrict C);  // A=B*C^{dag}
void times_dag12_U1(U1 * restrict A,
                    U1 const * restrict B,
                    U1 const * restrict  C); // A=B^{dag}*C^{dag}

void rand_matrix_U1(U1 * restrict A);

double norm_U1(U1 const * restrict A);

double retr_U1(U1 const * restrict A);
double imtr_U1(U1 const * restrict A);

void unitarize_U1(U1 * restrict A);
void taexp_U1(U1 * restrict A);

void print_on_screen_U1(U1 const * const A);
void print_on_file_U1(FILE *fp, U1 const * const A);
void print_on_binary_file_noswap_U1(FILE *fp, U1 const * const A);
void print_on_binary_file_swap_U1(FILE *fp, U1 const * const A);
void print_on_binary_file_bigen_U1(FILE *fp, U1 const * const A);
void read_from_file_U1(FILE *fp, U1 *A);
void read_from_binary_file_noswap_U1(FILE *fp, U1 *A);
void read_from_binary_file_swap_U1(FILE *fp, U1 *A);
void read_from_binary_file_bigen_U1(FILE *fp, U1 *A);

#endif // U1_H

