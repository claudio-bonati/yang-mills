#ifndef SU2_H
#define SU2_H

#include<stdio.h>

#include"macro.h"
#include"tens_prod.h"

//
// An Su2 matrix is represented as comp[0]+i\sum_{j=1}^3 comp[j]\sigma_j where
// sigma_j are Pauli matrices, comp[j] are real and \sum_{j=0}^3 comp[j]^2=1
//

typedef struct Su2 {
     double comp[4] __attribute__((aligned(DOUBLE_ALIGN)));
} Su2;

void init_Su2(Su2 * restrict A, double vec[4]);
void one_Su2(Su2 * restrict A);         // A=1
void zero_Su2(Su2 * restrict A);        // A=0

void equal_Su2(Su2 * restrict A, Su2 const * restrict B);      // A=B
void equal_dag_Su2(Su2 * restrict A, Su2 const * restrict B);  // A=B^{dag}

void plus_equal_Su2(Su2 * restrict A, Su2 const * restrict B);     // A+=B
void plus_equal_dag_Su2(Su2 * restrict A, Su2 const * restrict B); // A+=B^{dag}

void minus_equal_Su2(Su2 * restrict A, Su2 const * restrict B);     // A-=B
void minus_equal_times_real_Su2(Su2 * restrict A, Su2 const * restrict B, double r);     // A-=(r*B)
void minus_equal_dag_Su2(Su2 * restrict A, Su2 const * restrict B); // A-=B^{dag}

void lin_comb_Su2(Su2 * restrict A,
                  double b, Su2 const * restrict B,
                  double c, Su2 const * restrict C);       // A=b*B+c*C
void lin_comb_dag1_Su2(Su2 * restrict A,
                       double b, Su2 const * restrict B,
                       double c, Su2 const * restrict C);  // A=b*B^{dag}+c*C
void lin_comb_dag2_Su2(Su2 * restrict A,
                       double b, Su2 const * restrict B,
                       double c, Su2 const * restrict C);  // A=b*B+c*C^{dag}
void lin_comb_dag12_Su2(Su2 *A,
                        double b, Su2 const * restrict B,
                        double c, Su2 const * restrict C); // A=b*B^{dag}+c*C^{dag}

void times_equal_real_Su2(Su2 * restrict A, double r); // A*=r

void times_equal_Su2(Su2 * restrict A, Su2 const * restrict B);     // A*=B
void times_equal_dag_Su2(Su2 * restrict A, Su2 const * restrict B); // A*=B^{dag}

void times_Su2(Su2 * restrict A,
               Su2 const * restrict B,
               Su2 const * restrict C);       // A=B*C
void times_dag1_Su2(Su2 * restrict A,
                    Su2 const * restrict B,
                    Su2 const * restrict C);  // A=B^{dag}*C
void times_dag2_Su2(Su2 * restrict A,
                    Su2 const * restrict B,
                    Su2 const * restrict C);  // A=B*C^{dag}
void times_dag12_Su2(Su2 * restrict A,
                     Su2 const * restrict B,
                     Su2 const * restrict  C); // A=B^{dag}*C^{dag}

void rand_matrix_Su2(Su2 * restrict A);
void rand_matrix_p0_Su2(double p0, Su2 * restrict A);

double sqrtdet_Su2(Su2 const * restrict A);
double norm_Su2(Su2 const * restrict A);

double retr_Su2(Su2 const * restrict A);
double imtr_Su2(Su2 const * restrict A);

void unitarize_Su2(Su2 * restrict A);
void ta_Su2(Su2 * restrict A);
void taexp_Su2(Su2 * restrict A);

void print_on_screen_Su2(Su2 const * const A);
void print_on_file_Su2(FILE *fp, Su2 const * const A);
void print_on_binary_file_noswap_Su2(FILE *fp, Su2 const * const A);
void print_on_binary_file_swap_Su2(FILE *fp, Su2 const * const A);
void print_on_binary_file_bigen_Su2(FILE *fp, Su2 const * const A);
void read_from_file_Su2(FILE *fp, Su2 *A);
void read_from_binary_file_noswap_Su2(FILE *fp, Su2 *A);
void read_from_binary_file_swap_Su2(FILE *fp, Su2 *A);
void read_from_binary_file_bigen_Su2(FILE *fp, Su2 *A);

void TensProd_init_Su2(TensProd * restrict TP, Su2 const * restrict A1, Su2 const * restrict A2);

#endif
