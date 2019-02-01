#ifndef FUNCTION_POINTERS_C
#define FUNCTION_POINTERS_C

#include"../include/function_pointers.h"
#include"../include/macro.h"
#include"../include/su2.h"
#include"../include/su2_upd.h"
#include"../include/sun.h"
#include"../include/sun_upd.h"
#include"../include/tens_prod.h"
#include"../include/tens_prod_adj.h"
#include"../include/u1.h"
#include"../include/u1_upd.h"

#if NCOLOR == 1

void (*one)(GAUGE_GROUP *A)  = &one_U1;
void (*zero)(GAUGE_GROUP *A) = &zero_U1;

void (*equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &equal_U1;
void (*equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &equal_dag_U1;

void (*plus_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &plus_equal_U1;
void (*plus_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &plus_equal_dag_U1;

void (*minus_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &minus_equal_U1;
void (*minus_equal_times_real)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, double r) = &minus_equal_times_real_U1;
void (*minus_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &minus_equal_dag_U1;

void (*lin_comb)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_U1;
void (*lin_comb_dag1)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_dag1_U1;
void (*lin_comb_dag2)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_dag2_U1;
void (*lin_comb_dag12)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_dag12_U1;

void (*times_equal_real)(GAUGE_GROUP *A, double r) = &times_equal_real_U1;
void (*times_equal_complex)(GAUGE_GROUP *A, double complex r) = &times_equal_complex_U1;
void (*times_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &times_equal_U1;
void (*times_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const *B) = &times_equal_dag_U1;

void (*times)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_U1;
void (*times_dag1)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_dag1_U1;
void (*times_dag2)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_dag2_U1;
void (*times_dag12)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_dag12_U1;

void (*rand_matrix)(GAUGE_GROUP *A) = &rand_matrix_U1;

double (*norm)(GAUGE_GROUP const * const A) = &norm_U1;

double (*retr)(GAUGE_GROUP const * const A) = &retr_U1;
double (*imtr)(GAUGE_GROUP const * const A) = &imtr_U1;

void (*unitarize)(GAUGE_GROUP *A) = &unitarize_U1;
void (*ta)(GAUGE_GROUP *A) = &ta_U1;
void (*taexp)(GAUGE_GROUP *A) = &taexp_U1;

void (*print_on_screen)(GAUGE_GROUP const * const A) = &print_on_screen_U1;
void (*print_on_file)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_file_U1;
void (*print_on_binary_file_bigen)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_binary_file_bigen_U1;
void (*read_from_file)(FILE *fp, GAUGE_GROUP *A) = &read_from_file_U1;
void (*read_from_binary_file_bigen)(FILE *fp, GAUGE_GROUP *A) = &read_from_binary_file_bigen_U1;

void (*TensProd_init)(TensProd *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProd_init_U1;
void (*TensProdAdj_init)(TensProdAdj *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProdAdj_init_U1;

void (*single_heatbath)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple, GParam const * const param) = &single_heatbath_U1;
void (*single_overrelaxation)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &single_overrelaxation_U1;
void (*cool)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &cool_U1;

#elif NCOLOR == 2

void (*one)(GAUGE_GROUP *A)  = &one_Su2;
void (*zero)(GAUGE_GROUP *A) = &zero_Su2;

void (*equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &equal_Su2;
void (*equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &equal_dag_Su2;

void (*plus_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &plus_equal_Su2;
void (*plus_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &plus_equal_dag_Su2;

void (*minus_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &minus_equal_Su2;
void (*minus_equal_times_real)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, double r) = &minus_equal_times_real_Su2;
void (*minus_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &minus_equal_dag_Su2;

void (*lin_comb)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_Su2;
void (*lin_comb_dag1)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_dag1_Su2;
void (*lin_comb_dag2)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_dag2_Su2;
void (*lin_comb_dag12)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_dag12_Su2;

void (*times_equal_real)(GAUGE_GROUP *A, double r) = &times_equal_real_Su2;
void (*times_equal_complex)(GAUGE_GROUP *A, double complex r) = &times_equal_complex_Su2;
void (*times_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &times_equal_Su2;
void (*times_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const *B) = &times_equal_dag_Su2;

void (*times)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_Su2;
void (*times_dag1)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_dag1_Su2;
void (*times_dag2)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_dag2_Su2;
void (*times_dag12)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_dag12_Su2;

void (*rand_matrix)(GAUGE_GROUP *A) = &rand_matrix_Su2;

double (*norm)(GAUGE_GROUP const * const A) = &norm_Su2;

double (*retr)(GAUGE_GROUP const * const A) = &retr_Su2;
double (*imtr)(GAUGE_GROUP const * const A) = &imtr_Su2;

void (*unitarize)(GAUGE_GROUP *A) = &unitarize_Su2;
void (*ta)(GAUGE_GROUP *A) = &ta_Su2;
void (*taexp)(GAUGE_GROUP *A) = &taexp_Su2;

void (*print_on_screen)(GAUGE_GROUP const * const A) = &print_on_screen_Su2;
void (*print_on_file)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_file_Su2;
void (*print_on_binary_file_bigen)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_binary_file_bigen_Su2;
void (*read_from_file)(FILE *fp, GAUGE_GROUP *A) = &read_from_file_Su2;
void (*read_from_binary_file_bigen)(FILE *fp, GAUGE_GROUP *A) = &read_from_binary_file_bigen_Su2;

void (*TensProd_init)(TensProd *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProd_init_Su2;
void (*TensProdAdj_init)(TensProdAdj *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProdAdj_init_Su2;

void (*single_heatbath)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple, GParam const * const param) = &single_heatbath_Su2;
void (*single_overrelaxation)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &single_overrelaxation_Su2;
void (*cool)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &cool_Su2;

#else

void (*one)(GAUGE_GROUP *A)  = &one_SuN;
void (*zero)(GAUGE_GROUP *A) = &zero_SuN;

void (*equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &equal_SuN;
void (*equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &equal_dag_SuN;

void (*plus_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &plus_equal_SuN;
void (*plus_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &plus_equal_dag_SuN;

void (*minus_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &minus_equal_SuN;
void (*minus_equal_times_real)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, double r) = &minus_equal_times_real_SuN;
void (*minus_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &minus_equal_dag_SuN;

void (*lin_comb)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_SuN;
void (*lin_comb_dag1)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_dag1_SuN;
void (*lin_comb_dag2)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_dag2_SuN;
void (*lin_comb_dag12)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_dag12_SuN;

void (*times_equal_real)(GAUGE_GROUP *A, double r) = &times_equal_real_SuN;
void (*times_equal_complex)(GAUGE_GROUP *A, double complex r) = &times_equal_complex_SuN;
void (*times_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &times_equal_SuN;
void (*times_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const *B) = &times_equal_dag_SuN;

void (*times)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_SuN;
void (*times_dag1)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_dag1_SuN;
void (*times_dag2)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_dag2_SuN;
void (*times_dag12)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_dag12_SuN;

void (*rand_matrix)(GAUGE_GROUP *A) = &rand_matrix_SuN;

double (*norm)(GAUGE_GROUP const * const A) = &norm_SuN;

double (*retr)(GAUGE_GROUP const * const A) = &retr_SuN;
double (*imtr)(GAUGE_GROUP const * const A) = &imtr_SuN;

void (*unitarize)(GAUGE_GROUP *A) = &unitarize_SuN;
void (*ta)(GAUGE_GROUP *A) = &ta_SuN;
void (*taexp)(GAUGE_GROUP *A) = &taexp_SuN;

void (*print_on_screen)(GAUGE_GROUP const * const A) = &print_on_screen_SuN;
void (*print_on_file)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_file_SuN;
void (*print_on_binary_file_bigen)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_binary_file_bigen_SuN;
void (*read_from_file)(FILE *fp, GAUGE_GROUP *A) = &read_from_file_SuN;
void (*read_from_binary_file_bigen)(FILE *fp, GAUGE_GROUP *A) = &read_from_binary_file_bigen_SuN;

void (*TensProd_init)(TensProd *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProd_init_SuN;
void (*TensProdAdj_init)(TensProdAdj *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProdAdj_init_SuN;

void (*single_heatbath)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple, GParam const * const param) = &single_heatbath_SuN;
void (*single_overrelaxation)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &single_overrelaxation_SuN;
void (*cool)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &cool_SuN;

#endif



#endif
