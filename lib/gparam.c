#ifndef GPARAM_C
#define GPARAM_C

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/endianness.h"
#include"../include/gparam.h"
#include"../include/macro.h"


// remove from input file white/empty lines and comments
// comments start with the carachter #
void remove_white_line_and_comments(FILE *input)
  {
  int temp_i;

  temp_i=getc(input);
  if(temp_i=='\n' || temp_i==' ' || temp_i=='\043') // scan for white lines and comments
    {
    ungetc(temp_i, input);

    temp_i=getc(input);
    if(temp_i=='\n' || temp_i==' ') // white line
      {
      do
       {
       temp_i=getc(input);
       }
      while(temp_i=='\n' || temp_i==' ');
      }
    ungetc(temp_i, input);

    temp_i=getc(input);
    if(temp_i=='\043')  // comment, \043 = ascii oct for #
      {
      do
       {
       temp_i=getc(input);
       }
      while(temp_i!='\n');
      }
    else
      {
      ungetc(temp_i, input);
      }

    remove_white_line_and_comments(input);
    }
  else
    {
    ungetc(temp_i, input);
    }
  }


void readinput(char *in_file, GParam *param)
    {
    FILE *input;
    char str[STD_STRING_LENGTH], temp_str[STD_STRING_LENGTH];
    double temp_d;
    int temp_i, i;
    int err, end=1;
    unsigned int temp_ui;

    // this is to avoid unnecessary checks in case the multilevel is not used
    for(i=0; i<NLEVELS; i++)
       {
       param->d_ml_step[i]=0;
       }

    // just to avoid possible mistakes with uninitialized stuff
    for(i=0; i<NCOLOR; i++)
       {
       param->d_h[i]=0.0;
       }
    param->d_theta=0.0;
    param->d_mon_meas=0; // if =1 monopole measures are performed
    param->d_higgs_beta=0.0;

    input=fopen(in_file, "r");  // open the input file
    if(input==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    else
      {
      while(end==1)   // slide the file
           {
           remove_white_line_and_comments(input);

           err=fscanf(input, "%s", str);
           if(err!=1)
             {
             fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
             printf("err=%d\n", err);
             exit(EXIT_FAILURE);
             }

           if(strncmp(str, "size", 4)==0)
             {
             for(i=0; i<STDIM; i++)
                {
                err=fscanf(input, "%d", &temp_i);
                if(err!=1)
                  {
                  fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                  exit(EXIT_FAILURE);
                  }
                param->d_sizeg[i]=temp_i;
                }
             }

           else if(strncmp(str, "beta", 4)==0)
                  { 
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_beta=temp_d;
                  }
           else if(strncmp(str, "htracedef", 9)==0)
                  {
                  int halfncolor;
                  if(NCOLOR==1)
                    {
                    halfncolor=1;
                    }
                  else
                    {
                    halfncolor=(int)floor(NCOLOR/2.0);
                    }

                  for(i=0; i<halfncolor; i++)
                     {
                     err=fscanf(input, "%lf", &temp_d);
                     if(err!=1)
                       {
                       fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                       exit(EXIT_FAILURE);
                       }
                     param->d_h[i]=temp_d;
                     }
                  }
           else if(strncmp(str, "theta", 5)==0)
                  {
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_theta=temp_d;
                  }
           else if(strncmp(str, "higgs_beta", 10)==0)
                  {
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_higgs_beta=temp_d;
                  }

           else if(strncmp(str, "sample", 6)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_sample=temp_i;
                  }
           else if(strncmp(str, "thermal", 7)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_thermal=temp_i;
                  }
           else if(strncmp(str, "overrelax", 9)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_overrelax=temp_i;
                  }
           else if(strncmp(str, "measevery", 9)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_measevery=temp_i;
                  }
          
           else if(strncmp(str, "monomeas", 8)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_mon_meas=temp_i;
                  }


           else if(strncmp(str, "start", 5)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_start=temp_i;
                  }
           else if(strncmp(str, "saveconf_back_every", 19)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_saveconf_back_every=temp_i;
                  }
           else if(strncmp(str, "saveconf_analysis_every", 23)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_saveconf_analysis_every=temp_i;
                  }

           else if(strncmp(str, "epsilon_metro", 13)==0)
                  {
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_epsilon_metro=temp_d;
                  }

           else if(strncmp(str, "coolsteps", 9)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_coolsteps=temp_i;
                  }
           else if(strncmp(str, "coolrepeat", 10)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_coolrepeat=temp_i;
                  }

           else if(strncmp(str, "gfstep", 6)==0)
                  {
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_gfstep=temp_d;
                  }

           else if(strncmp(str, "multihit", 8)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_multihit=temp_i;
                  }
           else if(strncmp(str, "ml_step", 7)==0)
                  {
                  for(i=0; i<NLEVELS; i++)
                     {
                     err=fscanf(input, "%d", &temp_i);
                     if(err!=1)
                       {
                       fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                       exit(EXIT_FAILURE);
                       }
                     param->d_ml_step[i]=temp_i;
                     }
                  }
           else if(strncmp(str, "ml_upd", 6)==0)
                  {
                  for(i=0; i<NLEVELS; i++)
                     {
                     err=fscanf(input, "%d", &temp_i);
                     if(err!=1)
                       {
                       fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                       exit(EXIT_FAILURE);
                       }
                     param->d_ml_upd[i]=temp_i;
                     }
                  }
           else if(strncmp(str, "ml_level0_repeat", 16)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_ml_level0_repeat=temp_i;
                  }
           else if(strncmp(str, "dist_poly", 9)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_dist_poly=temp_i;
                  }
           else if(strncmp(str, "transv_dist", 11)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_trasv_dist=temp_i;
                  }
           else if(strncmp(str, "plaq_dir", 8)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_plaq_dir[0]=temp_i;

                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_plaq_dir[1]=temp_i;
                  }

           else if(strncmp(str, "conf_file", 9)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_conf_file, temp_str);
                  }
           else if(strncmp(str, "higgs_conf_file", 15)==0)
                  {
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_higgs_conf_file, temp_str);
                  }

           else if(strncmp(str, "data_file", 9)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_data_file, temp_str);
                  }
	   // MONOPOLES file
	   else if(strncmp(str, "mon_file", 8)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_mon_file, temp_str);
                  }
           else if(strncmp(str, "log_file", 8)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_log_file, temp_str);
                  }
           else if(strncmp(str, "ml_file", 7)==0)
                  {
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_ml_file, temp_str);
                  }

           else if(strncmp(str, "randseed", 8)==0)
                  { 
                  err=fscanf(input, "%u", &temp_ui);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_randseed=temp_ui;
                  }

           else
             {
             fprintf(stderr, "Error: unrecognized option %s in the file %s (%s, %d)\n", str, in_file, __FILE__, __LINE__);
             exit(EXIT_FAILURE);
             }

           remove_white_line_and_comments(input);

           // check if the read line is the last one
           temp_i=getc(input);
           if(temp_i==EOF)
             {
             end=0;
             }
           else
             {
             ungetc(temp_i, input);
             }
           }

      fclose(input);

      // VARIOUS CHECKS
      if(param->d_ml_step[0]!=0)
        {
        if(param->d_sizeg[0] % param->d_ml_step[0] || param->d_sizeg[0] < param->d_ml_step[0])
          {
          fprintf(stderr, "Error: size[0] has to be divisible by ml_step[0] and satisfy ml_step[0]<=size[0] (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        for(i=1; i<NLEVELS; i++)
           {
           if(param->d_ml_step[i-1] % param->d_ml_step[i] || param->d_ml_step[i-1] <= param->d_ml_step[i])
             {
             fprintf(stderr, "Error: ml_step[%d] has to be divisible by ml_step[%d] and larger than it (%s, %d)\n", i-1, i, __FILE__, __LINE__);
             exit(EXIT_FAILURE);
             }
           }
        if(param->d_ml_step[NLEVELS-1]==1)
          {
          fprintf(stderr, "Error: ml_step[%d] has to be larger than 1 (%s, %d)\n", NLEVELS-1, __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        }

      #ifdef OPENMP_MODE
      for(i=0; i<STDIM; i++)
         {
         temp_i = param->d_sizeg[i] % 2;
         if(temp_i!=0)
           {
           fprintf(stderr, "Error: size[%d] is not even.\n", i);
           fprintf(stderr, "When using OpenMP all the sides of the lattice have to be even! (%s, %d)\n", __FILE__, __LINE__);
           exit(EXIT_FAILURE);
           }
         }
      #endif

      err=0;
      for(i=0; i<STDIM; i++)
         {
         if(param->d_sizeg[i]==1)
           {
           err=1;
           }
         }
      if(err==1)
        {
        fprintf(stderr, "Error: all sizes has to be larger than 1: the totally reduced case is not implemented! (%s, %d)\n", __FILE__, __LINE__);
        }
      }
    }


// initialize data file
void init_data_file(FILE **dataf, GParam const * const param)
  {
  int i;

  if(param->d_start==2)
    {
    *dataf=fopen(param->d_data_file, "r");
    if(*dataf!=NULL) // file exists
      {
      fclose(*dataf);
      *dataf=fopen(param->d_data_file, "a");
      }
    else
      {
      *dataf=fopen(param->d_data_file, "w");
      fprintf(*dataf, "%d ", STDIM);
      for(i=0; i<STDIM; i++)
         {
         fprintf(*dataf, "%d ", param->d_sizeg[i]);
         }
      fprintf(*dataf, "\n");
      }
    }
  else
    {
    *dataf=fopen(param->d_data_file, "w");
    fprintf(*dataf, "%d ", STDIM);
    for(i=0; i<STDIM; i++)
       {
       fprintf(*dataf, "%d ", param->d_sizeg[i]);
       }
    fprintf(*dataf, "\n");
    }
  fflush(*dataf);
  }


// initialize monopoles file
void init_mon_file(FILE **monof, GParam const * const param)
  {
  int i;

  if(param->d_start==2)
    {
    *monof=fopen(param->d_mon_file, "r");
    if(*monof!=NULL) // file exists
      {
      fclose(*monof);
      *monof=fopen(param->d_mon_file, "a");
      }
    else
      {
      *monof=fopen(param->d_mon_file, "w");
      fprintf(*monof, "%d ", STDIM);
      for(i=0; i<STDIM; i++)
         {
         fprintf(*monof, "%d ", param->d_sizeg[i]);
         }
      fprintf(*monof, "\n");
      }
    }
  else
    {
    *monof=fopen(param->d_mon_file, "w");
    fprintf(*monof, "%d ", STDIM);
    for(i=0; i<STDIM; i++)
       {
       fprintf(*monof, "%d ", param->d_sizeg[i]);
       }
    fprintf(*monof, "\n");
    }
  fflush(*monof);
  }



// print simulation parameters
void print_parameters_local(GParam const * const param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+-----------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_local |\n");
    fprintf(fp, "+-----------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_sizeg[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_sizeg[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "monopoles: %d\n", param->d_mon_meas);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
    fprintf(fp, "\n");

    fprintf(fp, "coolsteps:      %d\n", param->d_coolsteps);
    fprintf(fp, "coolrepeat:     %d\n", param->d_coolrepeat);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


// print simulation parameters
void print_parameters_polycorr(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+--------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_polycorr |\n");
    fprintf(fp, "+--------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_sizeg[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_sizeg[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif

    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "\n");

    fprintf(fp, "multihit:   %d\n", param->d_multihit);
    fprintf(fp, "levels for multileves: %d\n", NLEVELS);
    fprintf(fp, "multilevel steps: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_step[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "updates for levels: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_upd[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "dist_poly:  %d\n", param->d_dist_poly);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


// print simulation parameters
void print_parameters_polycorr_higgs(GParam * param, time_t time_start, time_t time_end, double acc)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+--------------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_polycorr_higgs |\n");
    fprintf(fp, "+--------------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "number of higgs fields: %d\n", NHIGGS);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_sizeg[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_sizeg[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "higgs_beta: %.10lf\n", param->d_higgs_beta);

    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "\n");

    fprintf(fp, "epsilon_metro: %.10lf\n", param->d_epsilon_metro);
    fprintf(fp, "metropolis acceptance: %.10lf\n", acc);
    fprintf(fp, "\n");

    fprintf(fp, "multihit:   %d\n", param->d_multihit);
    fprintf(fp, "levels for multileves: %d\n", NLEVELS);
    fprintf(fp, "multilevel steps: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_step[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "updates for levels: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_upd[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "dist_poly:  %d\n", param->d_dist_poly);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }



// print simulation parameters
void print_parameters_polycorr_long(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+-------------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_polycorr_long |\n");
    fprintf(fp, "+-------------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_sizeg[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_sizeg[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "\n");

    fprintf(fp, "multihit:   %d\n", param->d_multihit);
    fprintf(fp, "levels for multileves: %d\n", NLEVELS);
    fprintf(fp, "multilevel steps: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_step[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "updates for levels: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_upd[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "level0_repeat:   %d\n", param->d_ml_level0_repeat);
    fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


// print simulation parameters
void print_parameters_polycorr_higgs_long(GParam * param, time_t time_start, time_t time_end, double acc)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+-------------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_polycorr_higgs_long |\n");
    fprintf(fp, "+-------------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "number of higgs fields: %d\n", NHIGGS);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_sizeg[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_sizeg[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "higgs_beta: %.10lf\n", param->d_higgs_beta);
    fprintf(fp, "\n");

    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "\n");

    fprintf(fp, "epsilon_metro: %.10lf\n", param->d_epsilon_metro);
    fprintf(fp, "metropolis acceptance (different from zero only when the whole lattice is updated!): %.10lf\n", acc);
    fprintf(fp, "\n");

    fprintf(fp, "multihit:   %d\n", param->d_multihit);
    fprintf(fp, "levels for multileves: %d\n", NLEVELS);
    fprintf(fp, "multilevel steps: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_step[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "updates for levels: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_upd[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "level0_repeat:   %d\n", param->d_ml_level0_repeat);
    fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


// print simulation parameters
void print_parameters_spectrum(GParam const * const param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+--------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_spectrum |\n");
    fprintf(fp, "+--------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_sizeg[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_sizeg[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }



// print simulation parameters
void print_parameters_t0(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+--------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_t0 |\n");
    fprintf(fp, "+--------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_sizeg[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_sizeg[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "gfstep:    %lf\n", param->d_gfstep);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }



// print simulation parameters for the tracedef case
void print_parameters_tracedef(GParam const * const param, time_t time_start, time_t time_end, double acc)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+--------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_tracedef |\n");
    fprintf(fp, "+--------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_sizeg[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_sizeg[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "h: %.10lf ", param->d_h[0]);
    for(i=1; i<(int) floor(NCOLOR/2.0); i++)
       {
       fprintf(fp, "%.10lf ", param->d_h[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "monopoles: %d\n", param->d_mon_meas);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
    fprintf(fp, "\n");

    fprintf(fp, "epsilon_metro: %.10lf\n", param->d_epsilon_metro);
    fprintf(fp, "metropolis acceptance: %.10lf\n", acc);
    fprintf(fp, "\n");

    fprintf(fp, "coolsteps:      %d\n", param->d_coolsteps);
    fprintf(fp, "coolrepeat:     %d\n", param->d_coolrepeat);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


// print simulation parameters
void print_parameters_tube_disc(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+---------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_tube_disc |\n");
    fprintf(fp, "+---------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_sizeg[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_sizeg[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "\n");

    fprintf(fp, "multihit:   %d\n", param->d_multihit);
    fprintf(fp, "levels for multileves: %d\n", NLEVELS);
    fprintf(fp, "multilevel steps: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_step[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "updates for levels: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_upd[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
    fprintf(fp, "transv_dist: %d\n", param->d_trasv_dist);
    fprintf(fp, "plaq_dir: %d %d\n", param->d_plaq_dir[0], param->d_plaq_dir[1]);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


// print simulation parameters
void print_parameters_tube_disc_long(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+--------------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_tube_disc_long |\n");
    fprintf(fp, "+--------------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_sizeg[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_sizeg[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "\n");

    fprintf(fp, "multihit:   %d\n", param->d_multihit);
    fprintf(fp, "levels for multileves: %d\n", NLEVELS);
    fprintf(fp, "multilevel steps: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_step[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "updates for levels: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_upd[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "level0_repeat:   %d\n", param->d_ml_level0_repeat);
    fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
    fprintf(fp, "transv_dist: %d\n", param->d_trasv_dist);
    fprintf(fp, "plaq_dir: %d %d\n", param->d_plaq_dir[0], param->d_plaq_dir[1]);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


// print simulation parameters
void print_parameters_tube_conn(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+---------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_tube_conn |\n");
    fprintf(fp, "+---------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_sizeg[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_sizeg[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "\n");

    fprintf(fp, "multihit:   %d\n", param->d_multihit);
    fprintf(fp, "levels for multileves: %d\n", NLEVELS);
    fprintf(fp, "multilevel steps: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_step[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "updates for levels: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_upd[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
    fprintf(fp, "transv_dist: %d\n", param->d_trasv_dist);
    fprintf(fp, "plaq_dir: %d %d\n", param->d_plaq_dir[0], param->d_plaq_dir[1]);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


// print simulation parameters
void print_parameters_tube_conn_long(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+--------------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_tube_conn_long |\n");
    fprintf(fp, "+--------------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_sizeg[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_sizeg[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "\n");

    fprintf(fp, "multihit:   %d\n", param->d_multihit);
    fprintf(fp, "levels for multileves: %d\n", NLEVELS);
    fprintf(fp, "multilevel steps: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_step[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "updates for levels: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_upd[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "level0_repeat:   %d\n", param->d_ml_level0_repeat);
    fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
    fprintf(fp, "transv_dist: %d\n", param->d_trasv_dist);
    fprintf(fp, "plaq_dir: %d %d\n", param->d_plaq_dir[0], param->d_plaq_dir[1]);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }

// print simulation parameters for the higgs case
void print_parameters_higgs(GParam const * const param, time_t time_start, time_t time_end, double acc)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+-----------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_higgs |\n");
    fprintf(fp, "+-----------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "number of higgs fields: %d\n", NHIGGS);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_sizeg[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_sizeg[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta:       %.10lf\n", param->d_beta);
    fprintf(fp, "higgs_beta: %.10lf ", param->d_higgs_beta);
    fprintf(fp, "\n\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "\n");

    fprintf(fp, "epsilon_metro: %.10lf\n", param->d_epsilon_metro);
    fprintf(fp, "metropolis acceptance: %.10lf\n", acc);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


#endif

