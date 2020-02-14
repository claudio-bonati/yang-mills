#ifndef GAUGE_CONF_H
#define GAUGE_CONF_H

#include"macro.h"

#include<complex.h>
#include<openssl/md5.h>
#include<stdio.h>

#include"flavour_matrix.h"
#include"gparam.h"
#include"geometry.h"
#include"su2.h"
#include"sun.h"
#include"tens_prod.h"
#include"tens_prod_adj.h"
#include"u1.h"


typedef struct Gauge_Conf {
  long update_index;

  GAUGE_GROUP **lattice;       // [volume] [STDIM]
  GAUGE_GROUP ***clover_array; // [volume] [STDIM] [STDIM]

  // for computing the polyakov loop correlator with multilevel
  TensProd ***ml_polycorr;   // [NLEVELS] [d_size[0]/d_ml_step[i]] [space_vol]
  GAUGE_GROUP **loc_poly;    // [d_size[0]/d_ml_step[NLEVELS-1]] [space_vol] auxilliary vector to be used in the multilevel

  // for computing the polyakov loop correlator in the adjoint rep. with multilevel
  TensProdAdj ***ml_polycorradj;   // [NLEVELS] [d_size[0]/d_ml_step[i]] [space_vol]
  GAUGE_GROUP_ADJ **loc_polyadj;   // [d_size[0]/d_ml_step[NLEVELS-1]] [space_vol] auxilliary vector to be used in the multilevel

  // for the disconnected correlator for string width
  TensProd **ml_polyplaq;        // [NLEVELS] [only slice 0] [space_vol]
  TensProdAdj **ml_polyplaqadj;  // [NLEVELS] [only slice 0] [space_vol]  for the adjoint case
  double complex *loc_plaq;      // [only slice 0] [space_vol] auxilliary vector to be used in the multilevel

  // for the connected correlator for string width
  TensProd **ml_polyplaqconn;   // [NLEVELS] [only slice 0] [space_vol]
  GAUGE_GROUP *loc_plaqconn;    // [only slice 0][space_vol] auxilliary vector to be used in the multilevel

  // for higgs field & co
  GAUGE_VECS *higgs;    // [volume]
  FMatrix *Qh;          // [volume]
  double complex *Dh;   // [volume]
 
  // for Abelian projection & co
  double ***diag_proj; // [volume] [STDIM] [NCOLOR]
  double **u1_subg;    // [volume] [STDIM]
  double **uflag;      // [volume] [STDIM] this is used to check if the link has already been considered in the searching of the wraps
 
  } Gauge_Conf;


// in gauge_conf_def.c
void init_gauge_conf(Gauge_Conf *GC,
                     GParam const * const param);
void read_gauge_conf(Gauge_Conf *GC,
                     GParam const * const param);
void free_gauge_conf(Gauge_Conf *GC,
                     GParam const * const param);
void write_conf_on_file_with_name(Gauge_Conf const * const GC,
                                  GParam const * const param,
                                  char const * const namefile);
void write_conf_on_file(Gauge_Conf const * const GC,
                        GParam const * const param);
void write_conf_on_file_back(Gauge_Conf const * const GC,
                             GParam const * const param);
void init_gauge_conf_from_gauge_conf(Gauge_Conf *GC,
                                     Gauge_Conf const * const GC2,
                                     GParam const * const param);
void compute_md5sum_conf(char *res,        // the lenght is 2*MD5_DIGEST_LENGTH
                         Gauge_Conf const * const GC,
                         GParam const * const param);

void alloc_polycorr_stuff(Gauge_Conf *GC,
                          GParam const * const param);
void free_polycorr_stuff(Gauge_Conf *GC,
                         GParam const * const param);
void write_polycorr_on_file(Gauge_Conf const * const GC,
                            GParam const * const param,
                            int iteration);
void read_polycorr_from_file(Gauge_Conf const * const GC,
                             GParam const * const param,
                             int *iteration);
void compute_md5sum_polycorr(char *res,        // the lenght is 2*MD5_DIGEST_LENGTH
                             Gauge_Conf const * const GC,
                             GParam const * const param);

void alloc_polycorradj(Gauge_Conf *GC,
                       GParam const * const param);
void free_polycorradj(Gauge_Conf *GC,
                      GParam const * const param);
void write_polycorradj_on_file(Gauge_Conf const * const GC,
                               GParam const * const param,
                               int iteration);
void read_polycorradj_from_file(Gauge_Conf const * const GC,
                                GParam const * const param,
                                int *iteration);
void compute_md5sum_polycorradj(char *res,        // the lenght is 2*MD5_DIGEST_LENGTH
                                Gauge_Conf const * const GC,
                                GParam const * const param);

void alloc_tube_disc_stuff(Gauge_Conf *GC,
                           GParam const * const param);
void free_tube_disc_stuff(Gauge_Conf *GC,
                          GParam const * const param);
void write_tube_disc_stuff_on_file(Gauge_Conf const * const GC,
                                   GParam const * const param,
                                   int iteration);
void read_tube_disc_stuff_from_file(Gauge_Conf const * const GC,
                                    GParam const * const param,
                                    int *iteration);
void compute_md5sum_tube_disc_stuff(char *res,     // the lenght is 2*MD5_DIGEST_LENGTH
                                    Gauge_Conf const * const GC,
                                    GParam const * const param);

void alloc_tubeadj_disc_stuff(Gauge_Conf *GC,
                              GParam const * const param);
void free_tubeadj_disc_stuff(Gauge_Conf *GC,
                             GParam const * const param);
void write_tubeadj_disc_stuff_on_file(Gauge_Conf const * const GC,
                                   GParam const * const param,
                                   int iteration);
void read_tubeadj_disc_stuff_from_file(Gauge_Conf const * const GC,
                                    GParam const * const param,
                                    int *iteration);
void compute_md5sum_tubeadj_disc_stuff(char *res,    // the lenght is 2*MD5_DIGEST_LENGTH
                                       Gauge_Conf const * const GC,
                                       GParam const * const param);

void alloc_tube_conn_stuff(Gauge_Conf *GC,
                           GParam const * const param);
void free_tube_conn_stuff(Gauge_Conf *GC,
                          GParam const * const param);
void write_tube_conn_stuff_on_file(Gauge_Conf const * const GC,
                                   GParam const * const param,
                                   int iteration);
void read_tube_conn_stuff_from_file(Gauge_Conf const * const GC,
                                    GParam const * const param,
                                    int *iteration);
void compute_md5sum_tube_conn_stuff(char *res,  // the lenght is 2*MD5_DIGEST_LENGTH
                                    Gauge_Conf const * const GC,
                                    GParam const * const param);

void alloc_clover_array(Gauge_Conf *GC,
                        GParam const * const param);
void end_clover_array(Gauge_Conf *GC,
                      GParam const * const param);

void init_higgs_conf(Gauge_Conf *GC,
                     GParam const * const param);
void read_higgs_conf(Gauge_Conf *GC,
                     GParam const * const param);
void free_higgs_conf(Gauge_Conf *GC);
void write_higgs_on_file_with_name(Gauge_Conf const * const GC,
                                   GParam const * const param,
                                   char const * const namefile);
void write_higgs_on_file(Gauge_Conf const * const GC,
                         GParam const * const param);
void write_higgs_on_file_back(Gauge_Conf const * const GC,
                              GParam const * const param);
void compute_md5sum_higgs(char *res,        // the lenght is 2*MD5_DIGEST_LENGTH
                          Gauge_Conf const * const GC,
                          GParam const * const param);

void alloc_diag_proj_stuff(Gauge_Conf *GC,
                           GParam const * const param);
void free_diag_proj_stuff(Gauge_Conf *GC,
                          GParam const * const param);


// in gauge_conf_meas.c
double plaquettep(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  long r,
                  int i,
                  int j);
double complex plaquettep_complex(Gauge_Conf const * const GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  long r,
                                  int i,
                                  int j);
void plaquettep_matrix(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r,
                       int i,
                       int j,
                       GAUGE_GROUP *matrix);
void clover(Gauge_Conf const * const GC,
            Geometry const * const geo,
            GParam const * const param,
            long r,
            int i,
            int j,
            GAUGE_GROUP *M);
void plaquette(Gauge_Conf const * const GC,
               Geometry const * const geo,
               GParam const * const param,
               double *plaqs,
               double *plaqt);
void plaquette_fundadj(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *plaqsf,
                       double *plaqtf,
                       double *plaqsa,
                       double *plaqta);
void clover_disc_energy(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        double *energy);
void polyakov(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              double *repoly,
              double *impoly);
void polyakov_adj(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  double *repoly,
                  double *impoly);
void polyakov_for_tracedef(Gauge_Conf const * const GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           double *repoly,
                           double *impoly);
double loc_topcharge(Gauge_Conf const * const GC,
                     Geometry const * const geo,
                     GParam const * const param,
                     long r);
double topcharge(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param);

void perform_measures_localobs(Gauge_Conf const * const GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep,
                               FILE *monofilep);

void perform_measures_localobs_with_tracedef(Gauge_Conf const * const GC,
                                             Geometry const * const geo,
                                             GParam const * const param,
                                             FILE *datafilep,
                                             FILE *monofilep);

void perform_measures_localobs_fundadj(Gauge_Conf const * const GC,
                                       Geometry const * const geo,
                                       GParam const * const param,
                                       FILE *datafilep);

void higgs_interaction(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *he);
void compute_flavour_observables(Gauge_Conf const * const GC,
                                 GParam const * const param,
                                 double *tildeG0,
                                 double *tildeGminp,
                                 double *tildeD0,
                                 double *tildeDminp);
void compute_flavour_observables_corr(Gauge_Conf const * const GC,
                                      Geometry const * const geo,
                                      GParam const * const param,
                                      double *corrQQ,
                                      double *corr0string0,
                                      double *corr0string1);
void perform_measures_higgs(Gauge_Conf * GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            FILE *datafilep);
void perform_measures_higgs_for_testing(Gauge_Conf *GC,
                                        Geometry const * const geo,
                                        GParam const * const param,
                                        FILE *datafilep);

void max_abelian_gauge_fix(Gauge_Conf *GC,
                           Geometry const * const geo,
                           GParam const * const param);
void diag_projection(Gauge_Conf *GC,
                     GParam const * const param);
void U1_extract(Gauge_Conf *GC,
                GParam const * const param,
                int subg);
void Di_Fjk(Gauge_Conf *GC,
            Geometry const * const geo,
            long r,
            int idir,
            int jdir,
            int kdir,
            double *DiFjk);
int DeGrand_current(Gauge_Conf *GC,
                    Geometry const * const geo,
                    long r,
                    int dir);
void wrap_search(Gauge_Conf *GC,
                 Geometry const * const geo,
                 GParam const * const param,
                 long r,
                 long r_tback,
                 int *num_wrap);
void monopoles_obs(Gauge_Conf *GC,
                   Geometry const * const geo,
                   GParam const * const param,
                   int subg,
                   FILE* monofilep);


// in gauge_conf_meas_multilevel.c
void optimize_multihit_polycorr(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep);
void optimize_multihit_polycorr_with_higgs(Gauge_Conf *GC,
                                           Geometry const * const geo,
                                           GParam const * const param,
                                           FILE *datafilep);
void optimize_multilevel_polycorr(Gauge_Conf *GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  FILE *datafilep);
void optimize_multilevel_polycorr_with_higgs(Gauge_Conf *GC,
                                             Geometry const * const geo,
                                             GParam const * const param,
                                             FILE *datafilep);
void perform_measures_polycorr(Gauge_Conf * GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep);
void perform_measures_polycorr_with_higgs(Gauge_Conf * GC,
                                          Geometry const * const geo,
                                          GParam const * const param,
                                          FILE *datafilep);

void optimize_multihit_polycorradj(Gauge_Conf *GC,
                                   Geometry const * const geo,
                                   GParam const * const param,
                                   FILE *datafilep);
void optimize_multilevel_polycorradj(Gauge_Conf *GC,
                                     Geometry const * const geo,
                                     GParam const * const param,
                                     FILE *datafilep);
void perform_measures_polycorradj(Gauge_Conf *GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  FILE *datafilep);

void optimize_multilevel_polycorr_long(Gauge_Conf *GC,
                                       GParam const * const param,
                                       FILE *datafilep);
void perform_measures_polycorr_long(Gauge_Conf * GC,
                                    GParam const * const param,
                                    FILE *datafilep);

void optimize_multilevel_polycorradj_long(Gauge_Conf *GC,
                                          GParam const * const param,
                                          FILE *datafilep);
void perform_measures_polycorradj_long(Gauge_Conf * GC,
                                       GParam const * const param,
                                       FILE *datafilep);

void perform_measures_tube_disc(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep);
void perform_measures_tube_disc_long(Gauge_Conf *GC,
                                     GParam const * const param,
                                     FILE *datafilep);

void perform_measures_tubeadj_disc(Gauge_Conf *GC,
                                   Geometry const * const geo,
                                   GParam const * const param,
                                   FILE *datafilep);
void perform_measures_tubeadj_disc_long(Gauge_Conf *GC,
                                        GParam const * const param,
                                        FILE *datafilep);

void perform_measures_tube_conn(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep);
void perform_measures_tube_conn_long(Gauge_Conf *GC,
                                     GParam const * const param,
                                     FILE *datafilep);


// in gauge_conf_multilevel.c
void multihit(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              long r,
              int dir,
              int num_hit,
              GAUGE_GROUP *G);
void multihit_with_higgs(Gauge_Conf const * const GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         long r,
                         int dir,
                         int num_hit,
                         GAUGE_GROUP *G);
void multihitadj(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param,
                 long r,
                 int dir,
                 int num_hit,
                 GAUGE_GROUP_ADJ *G);

void update_for_multilevel(Gauge_Conf * GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           int level);
void update_for_multilevel_with_higgs(Gauge_Conf * GC,
                                      Geometry const * const geo,
                                      GParam const * const param,
                                      int level);

void compute_local_poly(Gauge_Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param);
void compute_local_polyadj(Gauge_Conf *GC,
                           Geometry const * const geo,
                           GParam const * const param);
void multilevel_polycorr(Gauge_Conf *GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         int dt);
void multilevel_polycorr_with_higgs(Gauge_Conf *GC,
                                    Geometry const * const geo,
                                    GParam const * const param,
                                    int dt);
void multilevel_polycorradj(Gauge_Conf * GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            int dt);
void multilevel_polycorr_long(Gauge_Conf * GC,
                              Geometry const * const geo,
                              GParam const * const param,
                              int dt,
                              int iteration);
void multilevel_polycorr_long_with_higgs(Gauge_Conf * GC,
                                         Geometry const * const geo,
                                         GParam const * const param,
                                         int dt,
                                         int iteration);

void multilevel_polycorradj_long(Gauge_Conf * GC,
                                 Geometry const * const geo,
                                 GParam const * const param,
                                 int dt,
                                 int iteration);

void compute_local_poly_and_plaq(Gauge_Conf *GC,
                                 Geometry const * const geo,
                                 GParam const * const param);
void multilevel_tube_disc(Gauge_Conf *GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          int dt);
void multilevel_tube_disc_long(Gauge_Conf * GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               int dt,
                               int iteration);

void compute_local_polyadj_and_plaq(Gauge_Conf *GC,
                                    Geometry const * const geo,
                                    GParam const * const param);
void multilevel_tubeadj_disc(Gauge_Conf *GC,
                             Geometry const * const geo,
                             GParam const * const param,
                             int dt);
void multilevel_tubeadj_disc_long(Gauge_Conf * GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  int dt,
                                  int iteration);
   
void compute_local_poly_plaq_and_plaqconn(Gauge_Conf *GC,
                                          Geometry const * const geo,
                                          GParam const * const param);
void multilevel_tube_conn(Gauge_Conf * GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          int dt);
void multilevel_tube_conn_long(Gauge_Conf * GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               int dt,
                               int iteration);

// in gauge_conf_upd.c
void calcstaples_wilson(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const gparam,
                        long r,
                        int i,
                        GAUGE_GROUP *M);
void calcstaples_wilson_nosum(Gauge_Conf const * const GC,
                              Geometry const * const geo,
                              GParam const * const gparam,
                              long r,
                              int i,
                              GAUGE_GROUP *M);
void calcstaples_tracedef(Gauge_Conf const * const GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          long r,
                          int i,
                          GAUGE_GROUP * M);
void compute_clovers(Gauge_Conf const * const GC,
                     Geometry const * const geo,
                     GParam const * const param,
                     int direction);
void calcstaples_with_topo(Gauge_Conf const * const GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           long r,
                           int i,
                           GAUGE_GROUP *M);
void calcstaples_with_topo_nosum(Gauge_Conf const * const GC,
                                 Geometry const * const geo,
                                 GParam const * const param,
                                 long r,
                                 int i,
                                 GAUGE_GROUP *M);

void heatbath(Gauge_Conf * GC,
              Geometry const * const geo,
              GParam const * const param,
              long r,
              int i);
void overrelaxation(Gauge_Conf * GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    long r,
                    int i);
int metropolis(Gauge_Conf *GC,
               Geometry const * const geo,
               GParam const * const param,
               long r,
               int i,
               int numhits);
int metropolis_with_tracedef(Gauge_Conf *GC,
                             Geometry const * const geo,
                             GParam const * const param,
                             long r,
                             int i,
                             int numhits);
int metropolis_fundadj(Gauge_Conf *GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r,
                       int i,
                       int numhits);

void update(Gauge_Conf *GC,
            Geometry const * const geo,
            GParam const * const param);
void update_with_trace_def(Gauge_Conf *GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           double *acc);
void update_fundadj(Gauge_Conf *GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    double *acc);

void cooling(Gauge_Conf *GC,
             Geometry const * const geo,
             GParam const * const param,
             int n);
void gradflow_RKstep(Gauge_Conf *GC,
                     Gauge_Conf *helper1,
                     Gauge_Conf *helper2,
                     Geometry const * const geo,
                     GParam const *const param,
                     double dt);
void ape_smearing(Gauge_Conf *GC,
                  Geometry const * const geo,
                  GParam const *const param,
                  double alpha,
                  int n);

void heatbath_with_higgs(Gauge_Conf *GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         long r,
                         int i);
void overrelaxation_with_higgs(Gauge_Conf *GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               long r,
                               int i);
void calcstaples_for_higgs(Gauge_Conf *GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           long r,
                           GAUGE_VECS *staple);
void overrelaxation_for_higgs(Gauge_Conf *GC,
                              Geometry const * const geo,
                              GParam const * const param,
                              long r);
int metropolis_for_higgs(Gauge_Conf *GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         long r);
void update_with_higgs(Gauge_Conf * GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *acc);


#endif
