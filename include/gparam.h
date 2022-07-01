#ifndef GPARAM_H
#define GPARAM_H

#include<stdio.h>
#include<time.h>

#include"macro.h"

typedef struct GParam {
  // lattice dimensions
  int d_sizeg[STDIM];

  // simulation parameters
  double d_beta;
  double d_h[NCOLOR]; // parameters for the trace deformation
  double d_theta;
  double d_higgs_beta;

  // simulation details
  int d_sample;
  int d_thermal;
  int d_overrelax;
  int d_measevery;
  int d_mon_meas; // if 1 monopoles measures are performed.

  // initialization & saving
  int d_start;
  int d_saveconf_back_every;
  int d_saveconf_analysis_every;

  // for metropolis
  double d_epsilon_metro;

  // for cooling in measures
  int d_coolsteps;
  int d_coolrepeat;

  // for gradient-flow evolution
  double d_gfstep;

  // for multilevel
  int d_multihit;
  int d_ml_step[NLEVELS];
  int d_ml_upd[NLEVELS];
  int d_ml_level0_repeat;
  int d_dist_poly;
  int d_trasv_dist;
  int d_plaq_dir[2];

  // output file names
  char d_conf_file[STD_STRING_LENGTH];
  char d_higgs_conf_file[STD_STRING_LENGTH];
  char d_data_file[STD_STRING_LENGTH];
  char d_mon_file[STD_STRING_LENGTH];
  char d_log_file[STD_STRING_LENGTH];
  char d_ml_file[STD_STRING_LENGTH];

  // random seed
  unsigned int d_randseed;

} GParam;


void remove_white_line_and_comments(FILE *input);
void readinput(char *in_file, GParam *param);

void init_data_file(FILE **dataf, GParam const * const param);
void init_mon_file(FILE **monof, GParam const * const param);

void print_parameters_local(GParam const * const param, time_t time_start, time_t time_end);

void print_parameters_polycorr(GParam * param, time_t time_start, time_t time_end);
void print_parameters_polycorr_higgs(GParam * param, time_t time_start, time_t time_end, double acc);
void print_parameters_polycorr_long(GParam * param, time_t time_start, time_t time_end);
void print_parameters_polycorr_higgs_long(GParam * param, time_t time_start, time_t time_end, double acc);

void print_parameters_spectrum(GParam const * const param, time_t time_start, time_t time_end);

void print_parameters_t0(GParam * param, time_t time_start, time_t time_end);

void print_parameters_tracedef(GParam const * const param, time_t time_start, time_t time_end, double acc);

void print_parameters_tube_disc(GParam * param, time_t time_start, time_t time_end);
void print_parameters_tube_disc_long(GParam * param, time_t time_start, time_t time_end);

void print_parameters_tube_conn(GParam * param, time_t time_start, time_t time_end);
void print_parameters_tube_conn_long(GParam * param, time_t time_start, time_t time_end);

void print_parameters_higgs(GParam const * const param, time_t time_start, time_t time_end, double acc);

#endif
