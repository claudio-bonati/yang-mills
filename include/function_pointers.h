#ifndef FUNCTION_POINTERS_H
#define FUNCTION_POINTERS_H

#include<complex.h>
#include<stdio.h>

#include"macro.h"
#include"su2.h"
#include"su2_upd.h"
#include"sun.h"
#include"sun_upd.h"
#include"tens_prod.h"
#include"u1.h"
#include"u1_upd.h"

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
void (*print_on_file)(FILE *fp, GAUGE_GROUP const * const A);
void (*print_on_binary_file_bigen)(FILE *fp, GAUGE_GROUP const * const A);
void (*read_from_file)(FILE *fp, GAUGE_GROUP *A);
void (*read_from_binary_file_bigen)(FILE *fp, GAUGE_GROUP *A);
void (*read_from_binary_file_swap)(FILE *fp, GAUGE_GROUP *A);

void (*TensProd_init)(TensProd *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2);

void (*single_heatbath)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple, GParam const * const param);
void (*single_overrelaxation)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple);
void (*cool)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple);


// AUXILLIARY FUNCTIONS

void init_function_pointers(void);

#endif
