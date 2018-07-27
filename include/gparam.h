#ifndef GPARAM_H
#define GPARAM_H

#include<stdio.h>
#include<time.h>

#include"macro.h"

typedef struct GParam {
  // lattice dimensions
  int d_stdim; // initialized to STDIM in readinput;
  int d_size[STDIM];

  // simulation parameters
  double d_beta;
  double d_h[NCOLOR]; // parameters for the trace deformation

  // simulation details
  int d_sample;
  int d_thermal;
  int d_overrelax;
  int d_measevery;

  // initialization & saving
  int d_start;
  int d_saveconf_back_every;
  int d_saveconf_analysis_every;

  // for metropolis
  double d_epsilon_metro;

  //for 4d_from_5d
  double d_aniso_t;

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

  // output file names
  char d_conf_file[STD_STRING_LENGTH];
  char d_data_file[STD_STRING_LENGTH];
  char d_log_file[STD_STRING_LENGTH];
  char d_ml_file[STD_STRING_LENGTH];

  // random seed
  unsigned int d_randseed;

  // derived constants
  long d_volume;           // total volume
  double d_inv_vol;        // 1 / tot. volume
  long d_space_vol;        // spatial component of the volume
  double d_inv_space_vol;  // 1 / spatial component of the volume
} GParam;


void remove_white_line_and_comments(FILE *input);
void readinput(char *in_file, GParam *param);
void init_derived_constants(GParam *param);
void init_data_file(FILE **dataf, GParam const * const param);
void print_parameters_local(GParam const * const param, time_t time_start, time_t time_end);
void print_parameters_polycorr(GParam * param, time_t time_start, time_t time_end);
void print_parameters_polycorr_long(GParam * param, time_t time_start, time_t time_end);
void print_parameters_tube_disc(GParam * param, time_t time_start, time_t time_end);
void print_parameters_t0(GParam * param, time_t time_start, time_t time_end);
void print_parameters_tracedef(GParam const * const param, time_t time_start, time_t time_end, double acc);
void print_parameters_4d_from_5d(GParam const * const param, time_t time_start, time_t time_end);

#endif
