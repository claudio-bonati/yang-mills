#ifndef FUNCTION_POINTERS_H
#define FUNCTION_POINTERS_H

#include<complex.h>
#include<stdio.h>

#include"flavour_matrix.h"
#include"macro.h"
#include"son.h"
#include"son_upd.h"
#include"su2.h"
#include"su2_upd.h"
#include"sun.h"
#include"sun_upd.h"
#include"tens_prod.h"
#include"tens_prod_adj.h"
#include"u1.h"
#include"u1_upd.h"
#include"geometry.h"
#include"gparam.h"
#include"gauge_conf.h"

void (*one)(GAUGE_GROUP *A);         // A=1
void (*zero)(GAUGE_GROUP *A);        // A=0

void (*equal)(GAUGE_GROUP *A,
              GAUGE_GROUP const * const B);      // A=B
void (*equal_dag)(GAUGE_GROUP *A,
                  GAUGE_GROUP const * const B);  // A=B^{dag}

void (*plus_equal)(GAUGE_GROUP *A,
                   GAUGE_GROUP const * const B);     // A+=B
void (*plus_equal_dag)(GAUGE_GROUP *A,
                       GAUGE_GROUP const * const B); // A+=B^{dag}

void (*minus_equal)(GAUGE_GROUP *A,
                    GAUGE_GROUP const * const B);     // A-=B
void (*minus_equal_times_real)(GAUGE_GROUP *A,
                               GAUGE_GROUP const * const B, double r);   // A-=(r*B)
void (*minus_equal_dag)(GAUGE_GROUP *A,
                        GAUGE_GROUP const * const B); // A-=B^{dag}

void (*lin_comb)(GAUGE_GROUP *A,
                 double b, GAUGE_GROUP const * const B,
                 double c, GAUGE_GROUP const * const C);       // A=b*B+c*C
void (*lin_comb_dag1)(GAUGE_GROUP *A,
                      double b, GAUGE_GROUP const * const B,
                      double c, GAUGE_GROUP const * const C);  // A=b*B^{dag}+c*C
void (*lin_comb_dag2)(GAUGE_GROUP *A,
                      double b, GAUGE_GROUP const * const B,
                      double c, GAUGE_GROUP const * const C);  // A=b*B+c*C^{dag}
void (*lin_comb_dag12)(GAUGE_GROUP *A,
                       double b, GAUGE_GROUP const * const B,
                       double c, GAUGE_GROUP const * const C); // A=b*B^{dag}+c*C^{dag}

void (*times_equal_real)(GAUGE_GROUP *A, double r); // A*=r
void (*times_equal_complex)(GAUGE_GROUP *A, double complex r); // A*=r

void (*times_equal)(GAUGE_GROUP *A,
                    GAUGE_GROUP const * const B);     // A*=B
void (*times_equal_dag)(GAUGE_GROUP *A,
                        GAUGE_GROUP const *B); // A*=B^{dag}

void (*times)(GAUGE_GROUP *A,
              GAUGE_GROUP const * const B,
              GAUGE_GROUP const * const C);       // A=B*C
void (*times_dag1)(GAUGE_GROUP *A,
                   GAUGE_GROUP const * const B,
                   GAUGE_GROUP const * const C);  // A=B^{dag}*C
void (*times_dag2)(GAUGE_GROUP *A,
                   GAUGE_GROUP const * const B,
                   GAUGE_GROUP const * const C);  // A=B*C^{dag}
void (*times_dag12)(GAUGE_GROUP *A,
                    GAUGE_GROUP const * const B,
                    GAUGE_GROUP const * const C); // A=B^{dag}*C^{dag}

void (*rand_matrix)(GAUGE_GROUP *A);

double (*norm)(GAUGE_GROUP const * const A);

double (*retr)(GAUGE_GROUP const * const A);
double (*imtr)(GAUGE_GROUP const * const A);

void (*unitarize)(GAUGE_GROUP *A);
void (*ta)(GAUGE_GROUP *A);
void (*taexp)(GAUGE_GROUP *A);

void (*print_on_screen)(GAUGE_GROUP const * const A);
int  (*print_on_file)(FILE *fp, GAUGE_GROUP const * const A);
int  (*print_on_binary_file_bigen)(FILE *fp, GAUGE_GROUP const * const A);
int  (*read_from_file)(FILE *fp, GAUGE_GROUP *A);
int  (*read_from_binary_file_bigen)(FILE *fp, GAUGE_GROUP *A);

void (*diag_matrix_times)(GAUGE_GROUP * restrict A, double const lambda[NCOLOR], GAUGE_GROUP const * const restrict B);
void (*diag_matrix_times_dag)(GAUGE_GROUP * restrict A, double const lambda[NCOLOR], GAUGE_GROUP const * const restrict B);
void (*comp_MAG_gauge_transformation) (GAUGE_GROUP X_links[2*STDIM], double const lambda[NCOLOR], double OverRelaxParam, GAUGE_GROUP *G_mag);
void (*comp_outdiagnorm_of_X) (GAUGE_GROUP X_links[2*STDIM], double const lambda[NCOLOR], double *non_diag_contr);
void (*comp_functional_fmag) (GAUGE_GROUP X_links[2*STDIM], double const lambda[NCOLOR], double *fmag);
void (*diag_projection_single_site) (Gauge_Conf *GC, GAUGE_GROUP *link, long r, int dir);

void (*fund_to_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP const * const restrict B);

void (*TensProd_init)(TensProd *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2);
void (*TensProdAdj_init)(TensProdAdj *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2);
void (*TensProdAdj_initadj)(TensProdAdj *TP, GAUGE_GROUP_ADJ const * const A1, GAUGE_GROUP_ADJ const * const A2);

void (*one_adj)(GAUGE_GROUP_ADJ * restrict A);
void (*zero_adj)(GAUGE_GROUP_ADJ * restrict A);
void (*plus_equal_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP_ADJ const * const restrict B);
void (*times_equal_real_adj)(GAUGE_GROUP_ADJ * restrict A, double r);
void (*times_equal_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP_ADJ const * const restrict B);
double (*retr_adj)(GAUGE_GROUP_ADJ * restrict A);

void (*single_heatbath)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple);
void (*single_overrelaxation)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple);
void (*cool)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple);

void (*one_vecs)(GAUGE_VECS * restrict A);
void (*zero_vecs)(GAUGE_VECS * restrict A);

void (*equal_vecs)(GAUGE_VECS * restrict A,
                   GAUGE_VECS const * const restrict B);
void (*conjugate_vecs)(GAUGE_VECS * restrict A);
void (*minus_equal_vecs)(GAUGE_VECS * restrict A,
                         GAUGE_VECS const * const restrict B);
void (*plus_equal_vecs)(GAUGE_VECS * restrict A,
                        GAUGE_VECS const * const restrict B);
void (*times_equal_real_vecs)(GAUGE_VECS * restrict A, double r);
void (*times_equal_real_single_vecs)(GAUGE_VECS * restrict A, double r, int j);
void (*times_equal_complex_single_vecs)(GAUGE_VECS * restrict A, double complex r, int j);
double (*norm_vecs)(GAUGE_VECS const * const restrict A);
void (*normalize_vecs)(GAUGE_VECS * restrict A);

void (*rand_vecs)(GAUGE_VECS * restrict A);

double (*re_scal_prod_vecs)(GAUGE_VECS const * const restrict v1,
                            GAUGE_VECS const * const restrict v2);
double (*re_scal_prod_single_vecs)(GAUGE_VECS const * const restrict v1,
                                   GAUGE_VECS const * const restrict v2,
                                   int a,
                                   int b);
void (*matrix_times_vector_single_vecs)(GAUGE_VECS * restrict v1,
                                        GAUGE_GROUP const * const restrict matrix,
                                        GAUGE_VECS const * const restrict v2,
                                        int i);
void (*matrix_times_vector_all_vecs)(GAUGE_VECS * restrict v1,
                                     GAUGE_GROUP const * const restrict matrix,
                                     GAUGE_VECS const * const restrict v2);
void (*rotate_two_components_vecs)(GAUGE_VECS * restrict v1,
                                   GAUGE_VECS const * const restrict v2,
                                   int i,
                                   int j,
                                   double angle);
void (*vector_tensor_vector_vecs)(GAUGE_GROUP * restrict matrix,
                                  GAUGE_VECS const * const restrict v1,
                                  GAUGE_VECS const * const restrict v2);

void (*init_FMatrix_vecs)(FMatrix * restrict fmatrix,
                          GAUGE_VECS const * const restrict v1);
double complex (*HiggsU1Obs_vecs)(GAUGE_VECS const * const restrict v1);

int (*print_on_file_vecs)(FILE *fp, GAUGE_VECS const * const A);
int (*print_on_binary_file_bigen_vecs)(FILE *fp, GAUGE_VECS const * const A);
int (*read_from_file_vecs)(FILE *fp, GAUGE_VECS *A);
int (*read_from_binary_file_bigen_vecs)(FILE *fp, GAUGE_VECS *A);

void (*single_overrelaxation_vecs)(GAUGE_VECS *restrict link, GAUGE_VECS const * const staple);

#endif
