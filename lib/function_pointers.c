#ifndef FUNCTION_POINTERS_C
#define FUNCTION_POINTERS_C

#include"../include/macro.h"

#include"../include/function_pointers.h"

#include"../include/son.h"
#include"../include/son_monopoles.h"
#include"../include/son_upd.h"

#include"../include/su2.h"
#include"../include/su2_monopoles.h"
#include"../include/su2_upd.h"

#include"../include/sun.h"
#include"../include/sun_monopoles.h"
#include"../include/sun_upd.h"

#include"../include/tens_prod.h"
#include"../include/tens_prod_adj.h"

#include"../include/u1.h"
#include"../include/u1_monopoles.h"
#include"../include/u1_upd.h"

#if GGROUP == 0

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
int  (*print_on_file)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_file_U1;
int  (*print_on_binary_file_bigen)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_binary_file_bigen_U1;
int  (*read_from_file)(FILE *fp, GAUGE_GROUP *A) = &read_from_file_U1;
int  (*read_from_binary_file_bigen)(FILE *fp, GAUGE_GROUP *A) = &read_from_binary_file_bigen_U1;

void (*diag_matrix_times)(GAUGE_GROUP * restrict A, double const lambda[1], GAUGE_GROUP const * const restrict B) = &diag_matrix_times_U1;
void (*diag_matrix_times_dag)(GAUGE_GROUP * restrict A, double const lambda[1], GAUGE_GROUP const * const restrict B) = &diag_matrix_times_dag_U1;
void (*comp_outdiagnorm_of_X)(GAUGE_GROUP X_links[2*STDIM], double const lambda[1], double *non_diag_contr) = &comp_outdiagnorm_of_X_U1;
void (*comp_MAG_gauge_transformation) (GAUGE_GROUP X_links[2*STDIM], double const lambda[1], double OverRelaxParam, GAUGE_GROUP *G_mag) = &comp_MAG_gauge_transformation_U1;
void (*comp_functional_fmag) (GAUGE_GROUP X_links[2*STDIM], double const lambda[1], double *fmag) = &comp_functional_fmag_U1;
void (*diag_projection_single_site) (Gauge_Conf *GC, GAUGE_GROUP *link, long r, int dir) = &diag_projection_single_site_U1;

void (*fund_to_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP const * const restrict B)=&fund_to_adj_U1;

void (*TensProd_init)(TensProd *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProd_init_U1;
void (*TensProdAdj_init)(TensProdAdj *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProdAdj_init_U1;
void (*TensProdAdj_initadj)(TensProdAdj *TP, GAUGE_GROUP_ADJ const * const A1, GAUGE_GROUP_ADJ const * const A2) = &TensProdAdj_init_U1Adj;

void (*one_adj)(GAUGE_GROUP_ADJ * restrict A)=&one_U1Adj;
void (*zero_adj)(GAUGE_GROUP_ADJ * restrict A)=&zero_U1Adj;
void (*plus_equal_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP_ADJ const * const restrict B)=&plus_equal_U1Adj;
void (*times_equal_real_adj)(GAUGE_GROUP_ADJ * restrict A, double r)=&times_equal_real_U1Adj;
void (*times_equal_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP_ADJ const * const restrict B)=&times_equal_U1Adj;
double (*retr_adj)(GAUGE_GROUP_ADJ * restrict A)=&retr_U1Adj;

void (*single_heatbath)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &single_heatbath_U1;
void (*single_overrelaxation)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &single_overrelaxation_U1;
void (*cool)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &cool_U1;

void (*one_vecs)(GAUGE_VECS * restrict A)=&one_U1Vecs;
void (*zero_vecs)(GAUGE_VECS * restrict A)=&zero_U1Vecs;

void (*equal_vecs)(GAUGE_VECS * restrict A, GAUGE_VECS const * const restrict B)=&equal_U1Vecs;
void (*conjugate_vecs)(GAUGE_VECS * restrict A)=&conjugate_U1Vecs;
void (*minus_equal_vecs)(GAUGE_VECS * restrict A, GAUGE_VECS const * const restrict B)=&minus_equal_U1Vecs;
void (*plus_equal_vecs)(GAUGE_VECS * restrict A, GAUGE_VECS const * const restrict B)=&plus_equal_U1Vecs;
void (*times_equal_real_vecs)(GAUGE_VECS * restrict A, double r)=&times_equal_real_U1Vecs;
void (*times_equal_real_single_vecs)(GAUGE_VECS * restrict A, double r, int j)=&times_equal_real_single_U1Vecs;
void (*times_equal_complex_single_vecs)(GAUGE_VECS * restrict A, double complex r, int j)=&times_equal_complex_single_U1Vecs;

double (*norm_vecs)(GAUGE_VECS const * const restrict A)=&norm_U1Vecs;
void (*normalize_vecs)(GAUGE_VECS * restrict A)=&normalize_U1Vecs;

void (*rand_vecs)(GAUGE_VECS * restrict A)=&rand_vec_U1Vecs;

double (*re_scal_prod_vecs)(GAUGE_VECS const * const restrict v1, GAUGE_VECS const * const restrict v2)=&re_scal_prod_U1Vecs;
double (*re_scal_prod_single_vecs)(GAUGE_VECS const * const restrict v1, GAUGE_VECS const * const restrict v2, int a, int b)=&re_scal_prod_single_U1Vecs;
void (*matrix_times_vector_single_vecs)(GAUGE_VECS * restrict v1, GAUGE_GROUP const * const restrict matrix, GAUGE_VECS const * const restrict v2, int i)=&matrix_times_vector_single_U1Vecs;
void (*matrix_times_vector_all_vecs)(GAUGE_VECS * restrict v1, GAUGE_GROUP const * const restrict matrix, GAUGE_VECS const * const restrict v2)=&matrix_times_vector_all_U1Vecs;
void (*rotate_two_components_vecs)(GAUGE_VECS * restrict v1, GAUGE_VECS const * const restrict v2, int i, int j, double angle)=&rotate_two_components_U1Vecs;
void (*vector_tensor_vector_vecs)(GAUGE_GROUP * restrict matrix, GAUGE_VECS const * const restrict v1, GAUGE_VECS const * const restrict v2)=&vector_tensor_vector_U1Vecs;

void (*init_FMatrix_vecs)(FMatrix * restrict fmatrix, GAUGE_VECS const * const restrict v1)=&init_FMatrix_U1Vecs;
double complex (*HiggsU1Obs_vecs)(GAUGE_VECS const * const restrict v1)=&HiggsU1Obs_U1Vecs;

int (*print_on_file_vecs)(FILE *fp, GAUGE_VECS const * const A)=&print_on_file_U1Vecs;
int (*print_on_binary_file_bigen_vecs)(FILE *fp, GAUGE_VECS const * const A)=&print_on_binary_file_bigen_U1Vecs;
int (*read_from_file_vecs)(FILE *fp, GAUGE_VECS *A)=&read_from_file_U1Vecs;
int (*read_from_binary_file_bigen_vecs)(FILE *fp, GAUGE_VECS *A)=&read_from_binary_file_bigen_U1Vecs;

void (*single_overrelaxation_vecs)(GAUGE_VECS *restrict link, GAUGE_VECS const * const staple)=&single_overrelaxation_U1Vecs;

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
int  (*print_on_file)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_file_Su2;
int  (*print_on_binary_file_bigen)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_binary_file_bigen_Su2;
int  (*read_from_file)(FILE *fp, GAUGE_GROUP *A) = &read_from_file_Su2;
int  (*read_from_binary_file_bigen)(FILE *fp, GAUGE_GROUP *A) = &read_from_binary_file_bigen_Su2;

void (*diag_matrix_times)(GAUGE_GROUP * restrict A, double const lambda[2], GAUGE_GROUP const * const restrict B) = &diag_matrix_times_Su2;
void (*diag_matrix_times_dag)(GAUGE_GROUP * restrict A, double const lambda[2], GAUGE_GROUP const * const restrict B) = &diag_matrix_times_dag_Su2;
void (*comp_MAG_gauge_transformation) (GAUGE_GROUP helper_X[2*STDIM], double const lambda[2], double OverRelaxParam, GAUGE_GROUP *G_mag) = &comp_MAG_gauge_transformation_Su2;
void (*comp_outdiagnorm_of_X) (GAUGE_GROUP helper_X[2*STDIM], double const lambda[2], double *non_diag_contr) = &comp_outdiagnorm_of_X_Su2;
void (*comp_functional_fmag) (GAUGE_GROUP X_links[2*STDIM], double const lambda[2], double *fmag) = &comp_functional_fmag_Su2;
void (*diag_projection_single_site) (Gauge_Conf *GC, GAUGE_GROUP *link, long r, int dir) = &diag_projection_single_site_Su2;

void (*fund_to_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP const * const restrict B)=&fund_to_adj_Su2;

void (*TensProd_init)(TensProd *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProd_init_Su2;
void (*TensProdAdj_init)(TensProdAdj *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProdAdj_init_Su2;
void (*TensProdAdj_initadj)(TensProdAdj *TP, GAUGE_GROUP_ADJ const * const A1, GAUGE_GROUP_ADJ const * const A2) = &TensProdAdj_init_Su2Adj;

void (*one_adj)(GAUGE_GROUP_ADJ * restrict A)=&one_Su2Adj;
void (*zero_adj)(GAUGE_GROUP_ADJ * restrict A)=&zero_Su2Adj;
void (*plus_equal_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP_ADJ const * const restrict B)=&plus_equal_Su2Adj;
void (*times_equal_real_adj)(GAUGE_GROUP_ADJ * restrict A, double r)=&times_equal_real_Su2Adj;
void (*times_equal_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP_ADJ const * const restrict B)=&times_equal_Su2Adj;
double (*retr_adj)(GAUGE_GROUP_ADJ * restrict A)=&retr_Su2Adj;

void (*single_heatbath)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &single_heatbath_Su2;
void (*single_overrelaxation)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &single_overrelaxation_Su2;
void (*cool)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &cool_Su2;

void (*one_vecs)(GAUGE_VECS * restrict A)=&one_Su2Vecs;
void (*zero_vecs)(GAUGE_VECS * restrict A)=&zero_Su2Vecs;

void (*equal_vecs)(GAUGE_VECS * restrict A, GAUGE_VECS const * const restrict B)=&equal_Su2Vecs;
void (*conjugate_vecs)(GAUGE_VECS * restrict A)=&conjugate_Su2Vecs;
void (*minus_equal_vecs)(GAUGE_VECS * restrict A, GAUGE_VECS const * const restrict B)=&minus_equal_Su2Vecs;
void (*plus_equal_vecs)(GAUGE_VECS * restrict A, GAUGE_VECS const * const restrict B)=&plus_equal_Su2Vecs;
void (*times_equal_real_vecs)(GAUGE_VECS * restrict A, double r)=&times_equal_real_Su2Vecs;
void (*times_equal_real_single_vecs)(GAUGE_VECS * restrict A, double r, int j)=&times_equal_real_single_Su2Vecs;
void (*times_equal_complex_single_vecs)(GAUGE_VECS * restrict A, double complex r, int j)=&times_equal_complex_single_Su2Vecs;

double (*norm_vecs)(GAUGE_VECS const * const restrict A)=&norm_Su2Vecs;
void (*normalize_vecs)(GAUGE_VECS * restrict A)=&normalize_Su2Vecs;

void (*rand_vecs)(GAUGE_VECS * restrict A)=&rand_vec_Su2Vecs;

double (*re_scal_prod_vecs)(GAUGE_VECS const * const restrict v1, GAUGE_VECS const * const restrict v2)=&re_scal_prod_Su2Vecs;
double (*re_scal_prod_single_vecs)(GAUGE_VECS const * const restrict v1, GAUGE_VECS const * const restrict v2, int a, int b)=&re_scal_prod_single_Su2Vecs;
void (*matrix_times_vector_single_vecs)(GAUGE_VECS * restrict v1, GAUGE_GROUP const * const restrict matrix, GAUGE_VECS const * const restrict v2, int i)=&matrix_times_vector_single_Su2Vecs;
void (*matrix_times_vector_all_vecs)(GAUGE_VECS * restrict v1, GAUGE_GROUP const * const restrict matrix, GAUGE_VECS const * const restrict v2)=&matrix_times_vector_all_Su2Vecs;
void (*rotate_two_components_vecs)(GAUGE_VECS * restrict v1, GAUGE_VECS const * const restrict v2, int i, int j, double angle)=&rotate_two_components_Su2Vecs;
void (*vector_tensor_vector_vecs)(GAUGE_GROUP * restrict matrix, GAUGE_VECS const * const restrict v1, GAUGE_VECS const * const restrict v2)=&vector_tensor_vector_Su2Vecs;

void (*init_FMatrix_vecs)(FMatrix * restrict fmatrix, GAUGE_VECS const * const restrict v1)=&init_FMatrix_Su2Vecs;
double complex (*HiggsU1Obs_vecs)(GAUGE_VECS const * const restrict v1)=&HiggsU1Obs_Su2Vecs;

int (*print_on_file_vecs)(FILE *fp, GAUGE_VECS const * const A)=&print_on_file_Su2Vecs;
int (*print_on_binary_file_bigen_vecs)(FILE *fp, GAUGE_VECS const * const A)=&print_on_binary_file_bigen_Su2Vecs;
int (*read_from_file_vecs)(FILE *fp, GAUGE_VECS *A)=&read_from_file_Su2Vecs;
int (*read_from_binary_file_bigen_vecs)(FILE *fp, GAUGE_VECS *A)=&read_from_binary_file_bigen_Su2Vecs;

void (*single_overrelaxation_vecs)(GAUGE_VECS *restrict link, GAUGE_VECS const * const staple)=&single_overrelaxation_Su2Vecs;

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
int  (*print_on_file)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_file_SuN;
int  (*print_on_binary_file_bigen)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_binary_file_bigen_SuN;
int  (*read_from_file)(FILE *fp, GAUGE_GROUP *A) = &read_from_file_SuN;
int  (*read_from_binary_file_bigen)(FILE *fp, GAUGE_GROUP *A) = &read_from_binary_file_bigen_SuN;

void (*diag_matrix_times)(GAUGE_GROUP * restrict A, double const lambda[NCOLOR], GAUGE_GROUP const * const restrict B) = &diag_matrix_times_SuN;
void (*diag_matrix_times_dag)(GAUGE_GROUP * restrict A, double const lambda[NCOLOR], GAUGE_GROUP const * const restrict B) = &diag_matrix_times_dag_SuN;
void (*comp_outdiagnorm_of_X)(GAUGE_GROUP X_links[2*STDIM], double const lambda[NCOLOR], double *non_diag_contr) = &comp_outdiagnorm_of_X_SuN;
void (*comp_MAG_gauge_transformation) (GAUGE_GROUP X_links[2*STDIM], double const lambda[NCOLOR], double OverRelaxParam, GAUGE_GROUP *G_mag) = &comp_MAG_gauge_transformation_SuN;
void (*comp_functional_fmag) (GAUGE_GROUP X_links[2*STDIM], double const lambda[NCOLOR], double *fmag) = &comp_functional_fmag_SuN;
void (*diag_projection_single_site) (Gauge_Conf *GC, GAUGE_GROUP *link, long r, int dir) = &diag_projection_single_site_SuN;

void (*fund_to_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP const * const restrict B)=&fund_to_adj_SuN;

void (*TensProd_init)(TensProd *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProd_init_SuN;
void (*TensProdAdj_init)(TensProdAdj *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProdAdj_init_SuN;
void (*TensProdAdj_initadj)(TensProdAdj *TP, GAUGE_GROUP_ADJ const * const A1, GAUGE_GROUP_ADJ const * const A2) = &TensProdAdj_init_SuNAdj;

void (*one_adj)(GAUGE_GROUP_ADJ * restrict A)=&one_SuNAdj;
void (*zero_adj)(GAUGE_GROUP_ADJ * restrict A)=&zero_SuNAdj;
void (*plus_equal_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP_ADJ const * const restrict B)=&plus_equal_SuNAdj;
void (*times_equal_real_adj)(GAUGE_GROUP_ADJ * restrict A, double r)=&times_equal_real_SuNAdj;
void (*times_equal_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP_ADJ const * const restrict B)=&times_equal_SuNAdj;
double (*retr_adj)(GAUGE_GROUP_ADJ * restrict A)=&retr_SuNAdj;

void (*single_heatbath)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &single_heatbath_SuN;
void (*single_overrelaxation)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &single_overrelaxation_SuN;
void (*cool)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &cool_SuN;

void (*one_vecs)(GAUGE_VECS * restrict A)=&one_SuNVecs;
void (*zero_vecs)(GAUGE_VECS * restrict A)=&zero_SuNVecs;

void (*equal_vecs)(GAUGE_VECS * restrict A, GAUGE_VECS const * const restrict B)=&equal_SuNVecs;
void (*conjugate_vecs)(GAUGE_VECS * restrict A)=&conjugate_SuNVecs;
void (*minus_equal_vecs)(GAUGE_VECS * restrict A, GAUGE_VECS const * const restrict B)=&minus_equal_SuNVecs;
void (*plus_equal_vecs)(GAUGE_VECS * restrict A, GAUGE_VECS const * const restrict B)=&plus_equal_SuNVecs;
void (*times_equal_real_vecs)(GAUGE_VECS * restrict A, double r)=&times_equal_real_SuNVecs;
void (*times_equal_real_single_vecs)(GAUGE_VECS * restrict A, double r, int j)=&times_equal_real_single_SuNVecs;
void (*times_equal_complex_single_vecs)(GAUGE_VECS * restrict A, double complex r, int j)=&times_equal_complex_single_SuNVecs;

double (*norm_vecs)(GAUGE_VECS const * const restrict A)=&norm_SuNVecs;
void (*normalize_vecs)(GAUGE_VECS * restrict A)=&normalize_SuNVecs;

void (*rand_vecs)(GAUGE_VECS * restrict A)=&rand_vec_SuNVecs;

double (*re_scal_prod_vecs)(GAUGE_VECS const * const restrict v1, GAUGE_VECS const * const restrict v2)=&re_scal_prod_SuNVecs;
double (*re_scal_prod_single_vecs)(GAUGE_VECS const * const restrict v1, GAUGE_VECS const * const restrict v2, int a, int b)=&re_scal_prod_single_SuNVecs;
void (*matrix_times_vector_single_vecs)(GAUGE_VECS * restrict v1, GAUGE_GROUP const * const restrict matrix, GAUGE_VECS const * const restrict v2, int i)=&matrix_times_vector_single_SuNVecs;
void (*matrix_times_vector_all_vecs)(GAUGE_VECS * restrict v1, GAUGE_GROUP const * const restrict matrix, GAUGE_VECS const * const restrict v2)=&matrix_times_vector_all_SuNVecs;
void (*rotate_two_components_vecs)(GAUGE_VECS * restrict v1, GAUGE_VECS const * const restrict v2, int i, int j, double angle)=&rotate_two_components_SuNVecs;
void (*vector_tensor_vector_vecs)(GAUGE_GROUP * restrict matrix, GAUGE_VECS const * const restrict v1, GAUGE_VECS const * const restrict v2)=&vector_tensor_vector_SuNVecs;

void (*init_FMatrix_vecs)(FMatrix * restrict fmatrix, GAUGE_VECS const * const restrict v1)=&init_FMatrix_SuNVecs;
double complex (*HiggsU1Obs_vecs)(GAUGE_VECS const * const restrict v1)=&HiggsU1Obs_SuNVecs;

int (*print_on_file_vecs)(FILE *fp, GAUGE_VECS const * const A)=&print_on_file_SuNVecs;
int (*print_on_binary_file_bigen_vecs)(FILE *fp, GAUGE_VECS const * const A)=&print_on_binary_file_bigen_SuNVecs;
int (*read_from_file_vecs)(FILE *fp, GAUGE_VECS *A)=&read_from_file_SuNVecs;
int (*read_from_binary_file_bigen_vecs)(FILE *fp, GAUGE_VECS *A)=&read_from_binary_file_bigen_SuNVecs;

void (*single_overrelaxation_vecs)(GAUGE_VECS *restrict link, GAUGE_VECS const * const staple)=&single_overrelaxation_SuNVecs;

#endif

#elif GGROUP == 1

void (*one)(GAUGE_GROUP *A)  = &one_SoN;
void (*zero)(GAUGE_GROUP *A) = &zero_SoN;

void (*equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &equal_SoN;
void (*equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &equal_dag_SoN;

void (*plus_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &plus_equal_SoN;
void (*plus_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &plus_equal_dag_SoN;

void (*minus_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &minus_equal_SoN;
void (*minus_equal_times_real)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, double r) = &minus_equal_times_real_SoN;
void (*minus_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &minus_equal_dag_SoN;

void (*lin_comb)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_SoN;
void (*lin_comb_dag1)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_dag1_SoN;
void (*lin_comb_dag2)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_dag2_SoN;
void (*lin_comb_dag12)(GAUGE_GROUP *A, double b, GAUGE_GROUP const * const B, double c, GAUGE_GROUP const * const C) = &lin_comb_dag12_SoN;

void (*times_equal_real)(GAUGE_GROUP *A, double r) = &times_equal_real_SoN;
void (*times_equal_complex)(GAUGE_GROUP *A, double complex r) = &times_equal_complex_SoN;
void (*times_equal)(GAUGE_GROUP *A, GAUGE_GROUP const * const B) = &times_equal_SoN;
void (*times_equal_dag)(GAUGE_GROUP *A, GAUGE_GROUP const *B) = &times_equal_dag_SoN;

void (*times)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_SoN;
void (*times_dag1)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_dag1_SoN;
void (*times_dag2)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_dag2_SoN;
void (*times_dag12)(GAUGE_GROUP *A, GAUGE_GROUP const * const B, GAUGE_GROUP const * const C) = &times_dag12_SoN;

void (*rand_matrix)(GAUGE_GROUP *A) = &rand_matrix_SoN;

double (*norm)(GAUGE_GROUP const * const A) = &norm_SoN;

double (*retr)(GAUGE_GROUP const * const A) = &retr_SoN;
double (*imtr)(GAUGE_GROUP const * const A) = &imtr_SoN;

void (*unitarize)(GAUGE_GROUP *A) = &unitarize_SoN;
void (*ta)(GAUGE_GROUP *A) = &ta_SoN;
void (*taexp)(GAUGE_GROUP *A) = &taexp_SoN;

void (*print_on_screen)(GAUGE_GROUP const * const A) = &print_on_screen_SoN;
int  (*print_on_file)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_file_SoN;
int  (*print_on_binary_file_bigen)(FILE *fp, GAUGE_GROUP const * const A) = &print_on_binary_file_bigen_SoN;
int  (*read_from_file)(FILE *fp, GAUGE_GROUP *A) = &read_from_file_SoN;
int  (*read_from_binary_file_bigen)(FILE *fp, GAUGE_GROUP *A) = &read_from_binary_file_bigen_SoN;

void (*diag_matrix_times)(GAUGE_GROUP * restrict A, double const lambda[NCOLOR], GAUGE_GROUP const * const restrict B) = &diag_matrix_times_SoN;
void (*diag_matrix_times_dag)(GAUGE_GROUP * restrict A, double const lambda[NCOLOR], GAUGE_GROUP const * const restrict B) = &diag_matrix_times_dag_SoN;
void (*comp_outdiagnorm_of_X)(GAUGE_GROUP X_links[2*STDIM], double const lambda[NCOLOR], double *non_diag_contr) = &comp_outdiagnorm_of_X_SoN;
void (*comp_MAG_gauge_transformation) (GAUGE_GROUP X_links[2*STDIM], double const lambda[NCOLOR], double OverRelaxParam, GAUGE_GROUP *G_mag) = &comp_MAG_gauge_transformation_SoN;
void (*comp_functional_fmag) (GAUGE_GROUP X_links[2*STDIM], double const lambda[NCOLOR], double *fmag) = &comp_functional_fmag_SoN;
void (*diag_projection_single_site) (Gauge_Conf *GC, GAUGE_GROUP *link, long r, int dir) = &diag_projection_single_site_SoN;

void (*fund_to_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP const * const restrict B)=&fund_to_adj_SoN;

void (*TensProd_init)(TensProd *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProd_init_SoN;
void (*TensProdAdj_init)(TensProdAdj *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2) = &TensProdAdj_init_SoN;
void (*TensProdAdj_initadj)(TensProdAdj *TP, GAUGE_GROUP_ADJ const * const A1, GAUGE_GROUP_ADJ const * const A2) = &TensProdAdj_init_SoNAdj;

void (*one_adj)(GAUGE_GROUP_ADJ * restrict A)=&one_SoNAdj;
void (*zero_adj)(GAUGE_GROUP_ADJ * restrict A)=&zero_SoNAdj;
void (*plus_equal_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP_ADJ const * const restrict B)=&plus_equal_SoNAdj;
void (*times_equal_real_adj)(GAUGE_GROUP_ADJ * restrict A, double r)=&times_equal_real_SoNAdj;
void (*times_equal_adj)(GAUGE_GROUP_ADJ * restrict A, GAUGE_GROUP_ADJ const * const restrict B)=&times_equal_SoNAdj;
double (*retr_adj)(GAUGE_GROUP_ADJ * restrict A)=&retr_SoNAdj;

void (*single_heatbath)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &single_heatbath_SoN;
void (*single_overrelaxation)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &single_overrelaxation_SoN;
void (*cool)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple) = &cool_SoN;

void (*one_vecs)(GAUGE_VECS * restrict A)=&one_SoNVecs;
void (*zero_vecs)(GAUGE_VECS * restrict A)=&zero_SoNVecs;

void (*equal_vecs)(GAUGE_VECS * restrict A, GAUGE_VECS const * const restrict B)=&equal_SoNVecs;
void (*conjugate_vecs)(GAUGE_VECS * restrict A)=&conjugate_SoNVecs;
void (*minus_equal_vecs)(GAUGE_VECS * restrict A, GAUGE_VECS const * const restrict B)=&minus_equal_SoNVecs;
void (*plus_equal_vecs)(GAUGE_VECS * restrict A, GAUGE_VECS const * const restrict B)=&plus_equal_SoNVecs;
void (*times_equal_real_vecs)(GAUGE_VECS * restrict A, double r)=&times_equal_real_SoNVecs;
void (*times_equal_real_single_vecs)(GAUGE_VECS * restrict A, double r, int j)=&times_equal_real_single_SoNVecs;
void (*times_equal_complex_single_vecs)(GAUGE_VECS * restrict A, double complex r, int j)=&times_equal_complex_single_SoNVecs;

double (*norm_vecs)(GAUGE_VECS const * const restrict A)=&norm_SoNVecs;
void (*normalize_vecs)(GAUGE_VECS * restrict A)=&normalize_SoNVecs;

void (*rand_vecs)(GAUGE_VECS * restrict A)=&rand_vec_SoNVecs;

double (*re_scal_prod_vecs)(GAUGE_VECS const * const restrict v1, GAUGE_VECS const * const restrict v2)=&re_scal_prod_SoNVecs;
double (*re_scal_prod_single_vecs)(GAUGE_VECS const * const restrict v1, GAUGE_VECS const * const restrict v2, int a, int b)=&re_scal_prod_single_SoNVecs;
void (*matrix_times_vector_single_vecs)(GAUGE_VECS * restrict v1, GAUGE_GROUP const * const restrict matrix, GAUGE_VECS const * const restrict v2, int i)=&matrix_times_vector_single_SoNVecs;
void (*matrix_times_vector_all_vecs)(GAUGE_VECS * restrict v1, GAUGE_GROUP const * const restrict matrix, GAUGE_VECS const * const restrict v2)=&matrix_times_vector_all_SoNVecs;
void (*rotate_two_components_vecs)(GAUGE_VECS * restrict v1, GAUGE_VECS const * const restrict v2, int i, int j, double angle)=&rotate_two_components_SoNVecs;
void (*vector_tensor_vector_vecs)(GAUGE_GROUP * restrict matrix, GAUGE_VECS const * const restrict v1, GAUGE_VECS const * const restrict v2)=&vector_tensor_vector_SoNVecs;

void (*init_FMatrix_vecs)(FMatrix * restrict fmatrix, GAUGE_VECS const * const restrict v1)=&init_FMatrix_SoNVecs;
double complex (*HiggsU1Obs_vecs)(GAUGE_VECS const * const restrict v1)=&HiggsU1Obs_SoNVecs;

int (*print_on_file_vecs)(FILE *fp, GAUGE_VECS const * const A)=&print_on_file_SoNVecs;
int (*print_on_binary_file_bigen_vecs)(FILE *fp, GAUGE_VECS const * const A)=&print_on_binary_file_bigen_SoNVecs;
int (*read_from_file_vecs)(FILE *fp, GAUGE_VECS *A)=&read_from_file_SoNVecs;
int (*read_from_binary_file_bigen_vecs)(FILE *fp, GAUGE_VECS *A)=&read_from_binary_file_bigen_SoNVecs;

void (*single_overrelaxation_vecs)(GAUGE_VECS *restrict link, GAUGE_VECS const * const staple)=&single_overrelaxation_SoNVecs;


#endif

#endif
