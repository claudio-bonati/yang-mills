#ifndef GAUGE_CONF_H
#define GAUGE_CONF_H

#include"macro.h"

#ifdef HASH_MODE
  #include<openssl/md5.h>
#endif
#include<stdio.h>

#include"gparam.h"
#include"geometry.h"


typedef struct Gauge_Conf {
  long update_index;

  GAUGE_GROUP **lattice;

  // for computing the polyakov loop correlator with multilevel
  TensProd **ml_polycorr_ris;
  TensProd **ml_polycorr_tmp;

  // for the disconnected correlator for string width
  TensProd ***ml_polyplaq_ris;
  TensProd ***ml_polyplaq_tmp;

  } Gauge_Conf;


// in gauge_conf_def.c
void init_gauge_conf(Gauge_Conf *GC,
                     GParam const * const param);
void read_gauge_conf(Gauge_Conf *GC,
                     GParam const * const param);
void end_gauge_conf(Gauge_Conf *GC,
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
void init_4d_gauge_conf_from_5d_gauge_conf(Gauge_Conf *GC4d,
                                           Gauge_Conf const * const GC5d,
                                           GParam const * const param5d,
                                           int t);
void compute_md5sum(char *res,        // the lenght is 2*MD5_DIGEST_LENGTH
                    Gauge_Conf const * const GC,
                    GParam const * const param);

void init_polycorr(Gauge_Conf *GC,
                   GParam const * const param);
void end_polycorr(Gauge_Conf *GC);
void write_polycorr_on_file(Gauge_Conf const * const GC,
                            GParam const * const param,
                            int tstart,
                            int iteration);
void read_polycorr_from_file(Gauge_Conf const * const GC,
                             GParam const * const param,
                             int *tstart,
                             int *iteration);
void compute_md5sum_polycorr(char *res,        // the lenght is 2*MD5_DIGEST_LENGTH
                             Gauge_Conf const * const GC,
                             GParam const * const param);

void init_polycorr_and_polyplaq(Gauge_Conf *GC,
                                GParam const * const param);
void end_polycorr_and_polyplaq(Gauge_Conf *GC,
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
void clover(Gauge_Conf const * const GC,
            Geometry const * const geo,
            GParam const * const param,
            long r,
            int j,
            int i,
            GAUGE_GROUP *M);
void plaquette(Gauge_Conf const * const GC,
               Geometry const * const geo,
               GParam const * const param,
               double *plaqs,
               double *plaqt);
void clover_disc_energy(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        double *energy);
void polyakov(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              double *repoly,
              double *impoly);
double topcharge(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param);
void perform_measures_localobs(Gauge_Conf const * const GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep);
void optimize_multihit_polycorr(Gauge_Conf *GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       FILE *datafilep);
void optimize_multilevel_potQbarQ(Gauge_Conf *GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  FILE *datafilep);
void perform_measures_pot_QbarQ(Gauge_Conf * GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep);
void perform_measures_pot_QbarQ_long(Gauge_Conf * GC,
                                     GParam const * const param,
                                     FILE *datafilep);
void optimize_multilevel_tube_disc(Gauge_Conf *GC,
                                     Geometry const * const geo,
                                     GParam const * const param,
                                     FILE *datafilep);
void perform_measures_tube_disc(Gauge_Conf *GC,
                                   Geometry const * const geo,
                                   GParam const * const param,
                                   FILE *datafilep);
void perform_measures_4d_from_5d(Gauge_Conf const * const GC,
                                 Geometry const * const geo,
                                 GParam const * const param5d,
                                 FILE *datafilep);


// in gauge_conf_multilevel.c
void multihit(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              long r,
              int dir,
              int num_hit,
              GAUGE_GROUP *G);
void compute_local_poly(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        int t_start,
                        int dt,
                        GAUGE_GROUP *loc_poly);
void slice_single_update(Gauge_Conf *GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         int t_start,
                         int dt);
void multilevel_pot_QbarQ(Gauge_Conf *GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          int t_start,
                          int dt);
void multilevel_pot_QbarQ_long(Gauge_Conf * GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               int t_start,
                               int dt,
                               int iteration);
void compute_plaq_on_slice1(Gauge_Conf const * const GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            double complex **plaq);
void multilevel_tube_disc_QbarQ(Gauge_Conf * GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                int t_start,
                                int dt);


// in gauge_conf_upd.c
void calcstaples_wilson(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const gparam,
                        long r,
                        int i,
                        GAUGE_GROUP *M);
void calcstaples_wilson_totred(Gauge_Conf const * const GC,
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
void calcstaples_wilson_with_aniso_t(Gauge_Conf const * const GC,
                                     Geometry const * const geo,
                                     GParam const * const gparam,
                                     long r,
                                     int i,
                                     GAUGE_GROUP *M);
void heatbath(Gauge_Conf * GC,
              Geometry const * const geo,
              GParam const * const param,
              long r,
              int i);
void heatbath_with_aniso_t(Gauge_Conf * GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           long r,
                           int i);
void overrelaxation(Gauge_Conf * GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    long r,
                    int i);
void overrelaxation_with_aniso_t(Gauge_Conf * GC,
                                 Geometry const * const geo,
                                 GParam const * const param,
                                 long r,
                                 int i);
int metropolis(Gauge_Conf *GC,
               Geometry const * const geo,
               GParam const * const param,
               long r,
               int i);
int metropolis_with_tracedef(Gauge_Conf *GC,
                             Geometry const * const geo,
                             GParam const * const param,
                             long r,
                             int i);
int metropolis_totred(Gauge_Conf *GC,
                      Geometry const * const geo,
                      GParam const * const param,
                      long r,
                      int i);
int metropolis_with_tracedef_totred(Gauge_Conf *GC,
                                    Geometry const * const geo,
                                    GParam const * const param,
                                    long r,
                                    int i);
void update(Gauge_Conf * GC,
            Geometry const * const geo,
            GParam const * const param);
void update_with_trace_def(Gauge_Conf * GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           double *acc);
void update_with_trace_def_nototred(Gauge_Conf * GC,
                                    Geometry const * const geo,
                                    GParam const * const param,
                                    double *acc);
void update_with_trace_def_totred(Gauge_Conf * GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  double *acc);
void update_with_aniso_t(Gauge_Conf * GC,
                         Geometry const * const geo,
                         GParam const * const param);
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

#endif
