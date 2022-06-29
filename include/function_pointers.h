#ifndef FUNCTION_POINTERS_H
#define FUNCTION_POINTERS_H

#include<complex.h>
#include<stdio.h>

#include"flavour_matrix.h"
#include"macro.h"
#include"su2_upd.h"
#include"sun.h"
#include"sun_upd.h"
#include"tens_prod.h"
#include"geometry.h"
#include"gparam.h"
#include"gauge_conf.h"

extern void (*one)(GAUGE_GROUP *A);         // A=1
extern void (*zero)(GAUGE_GROUP *A);        // A=0

extern void (*equal)(GAUGE_GROUP *A,
              GAUGE_GROUP const * const B);      // A=B
extern void (*equal_dag)(GAUGE_GROUP *A,
                  GAUGE_GROUP const * const B);  // A=B^{dag}

extern void (*plus_equal)(GAUGE_GROUP *A,
                   GAUGE_GROUP const * const B);     // A+=B
extern void (*plus_equal_dag)(GAUGE_GROUP *A,
                       GAUGE_GROUP const * const B); // A+=B^{dag}

extern void (*minus_equal)(GAUGE_GROUP *A,
                    GAUGE_GROUP const * const B);     // A-=B
extern void (*minus_equal_times_real)(GAUGE_GROUP *A,
                               GAUGE_GROUP const * const B, double r);   // A-=(r*B)
extern void (*minus_equal_dag)(GAUGE_GROUP *A,
                        GAUGE_GROUP const * const B); // A-=B^{dag}

extern void (*lin_comb)(GAUGE_GROUP *A,
                 double b, GAUGE_GROUP const * const B,
                 double c, GAUGE_GROUP const * const C);       // A=b*B+c*C
extern void (*lin_comb_dag1)(GAUGE_GROUP *A,
                      double b, GAUGE_GROUP const * const B,
                      double c, GAUGE_GROUP const * const C);  // A=b*B^{dag}+c*C
extern void (*lin_comb_dag2)(GAUGE_GROUP *A,
                      double b, GAUGE_GROUP const * const B,
                      double c, GAUGE_GROUP const * const C);  // A=b*B+c*C^{dag}
extern void (*lin_comb_dag12)(GAUGE_GROUP *A,
                       double b, GAUGE_GROUP const * const B,
                       double c, GAUGE_GROUP const * const C); // A=b*B^{dag}+c*C^{dag}

extern void (*times_equal_real)(GAUGE_GROUP *A, double r); // A*=r
extern void (*times_equal_complex)(GAUGE_GROUP *A, double complex r); // A*=r

extern void (*times_equal)(GAUGE_GROUP *A,
                    GAUGE_GROUP const * const B);     // A*=B
extern void (*times_equal_dag)(GAUGE_GROUP *A,
                        GAUGE_GROUP const *B); // A*=B^{dag}

extern void (*times)(GAUGE_GROUP *A,
              GAUGE_GROUP const * const B,
              GAUGE_GROUP const * const C);       // A=B*C
extern void (*times_dag1)(GAUGE_GROUP *A,
                   GAUGE_GROUP const * const B,
                   GAUGE_GROUP const * const C);  // A=B^{dag}*C
extern void (*times_dag2)(GAUGE_GROUP *A,
                   GAUGE_GROUP const * const B,
                   GAUGE_GROUP const * const C);  // A=B*C^{dag}
extern void (*times_dag12)(GAUGE_GROUP *A,
                    GAUGE_GROUP const * const B,
                    GAUGE_GROUP const * const C); // A=B^{dag}*C^{dag}

extern void (*rand_matrix)(GAUGE_GROUP *A);

extern double (*norm)(GAUGE_GROUP const * const A);

extern double (*retr)(GAUGE_GROUP const * const A);
extern double (*imtr)(GAUGE_GROUP const * const A);

extern void (*unitarize)(GAUGE_GROUP *A);
extern void (*ta)(GAUGE_GROUP *A);
extern void (*taexp)(GAUGE_GROUP *A);

extern void (*print_on_screen)(GAUGE_GROUP const * const A);
extern int  (*print_on_file)(FILE *fp, GAUGE_GROUP const * const A);
extern int  (*print_on_binary_file_bigen)(FILE *fp, GAUGE_GROUP const * const A);
extern int  (*read_from_file)(FILE *fp, GAUGE_GROUP *A);
extern int  (*read_from_binary_file_bigen)(FILE *fp, GAUGE_GROUP *A);

extern void (*diag_matrix_times)(GAUGE_GROUP * restrict A, double const lambda[NCOLOR], GAUGE_GROUP const * const restrict B);
extern void (*diag_matrix_times_dag)(GAUGE_GROUP * restrict A, double const lambda[NCOLOR], GAUGE_GROUP const * const restrict B);
extern void (*comp_MAG_gauge_transformation) (GAUGE_GROUP X_links[2*STDIM], double const lambda[NCOLOR], double OverRelaxParam, GAUGE_GROUP *G_mag);
extern void (*comp_outdiagnorm_of_X) (GAUGE_GROUP X_links[2*STDIM], double const lambda[NCOLOR], double *non_diag_contr);
extern void (*comp_functional_fmag) (GAUGE_GROUP X_links[2*STDIM], double const lambda[NCOLOR], double *fmag);
extern void (*diag_projection_single_site) (Gauge_Conf *GC, GAUGE_GROUP *link, long r, int dir);

extern void (*TensProd_init)(TensProd *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2);

extern void (*single_heatbath)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple);
extern void (*single_overrelaxation)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple);
extern void (*cool)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple);

extern void (*one_vecs)(GAUGE_VECS * restrict A);
extern void (*zero_vecs)(GAUGE_VECS * restrict A);

extern void (*equal_vecs)(GAUGE_VECS * restrict A,
                   GAUGE_VECS const * const restrict B);
extern void (*conjugate_vecs)(GAUGE_VECS * restrict A);
extern void (*minus_equal_vecs)(GAUGE_VECS * restrict A,
                         GAUGE_VECS const * const restrict B);
extern void (*plus_equal_vecs)(GAUGE_VECS * restrict A,
                        GAUGE_VECS const * const restrict B);
extern void (*times_equal_real_vecs)(GAUGE_VECS * restrict A, double r);
extern void (*times_equal_real_single_vecs)(GAUGE_VECS * restrict A, double r, int j);
extern void (*times_equal_complex_single_vecs)(GAUGE_VECS * restrict A, double complex r, int j);
extern double (*norm_vecs)(GAUGE_VECS const * const restrict A);
extern void (*normalize_vecs)(GAUGE_VECS * restrict A);

extern void (*rand_vecs)(GAUGE_VECS * restrict A);

extern double (*re_scal_prod_vecs)(GAUGE_VECS const * const restrict v1,
                            GAUGE_VECS const * const restrict v2);
extern double (*re_scal_prod_single_vecs)(GAUGE_VECS const * const restrict v1,
                                   GAUGE_VECS const * const restrict v2,
                                   int a,
                                   int b);
extern void (*matrix_times_vector_single_vecs)(GAUGE_VECS * restrict v1,
                                        GAUGE_GROUP const * const restrict matrix,
                                        GAUGE_VECS const * const restrict v2,
                                        int i);
extern void (*matrix_times_vector_all_vecs)(GAUGE_VECS * restrict v1,
                                     GAUGE_GROUP const * const restrict matrix,
                                     GAUGE_VECS const * const restrict v2);
extern void (*rotate_two_components_vecs)(GAUGE_VECS * restrict v1,
                                   GAUGE_VECS const * const restrict v2,
                                   int i,
                                   int j,
                                   double angle);
extern void (*vector_tensor_vector_vecs)(GAUGE_GROUP * restrict matrix,
                                  GAUGE_VECS const * const restrict v1,
                                  GAUGE_VECS const * const restrict v2);

extern void (*init_FMatrix_vecs)(FMatrix * restrict fmatrix,
                          GAUGE_VECS const * const restrict v1);
extern double complex (*HiggsU1Obs_vecs)(GAUGE_VECS const * const restrict v1);

extern int (*print_on_file_vecs)(FILE *fp, GAUGE_VECS const * const A);
extern int (*print_on_binary_file_bigen_vecs)(FILE *fp, GAUGE_VECS const * const A);
extern int (*read_from_file_vecs)(FILE *fp, GAUGE_VECS *A);
extern int (*read_from_binary_file_bigen_vecs)(FILE *fp, GAUGE_VECS *A);

extern void (*single_overrelaxation_vecs)(GAUGE_VECS *restrict link, GAUGE_VECS const * const staple);

#endif
