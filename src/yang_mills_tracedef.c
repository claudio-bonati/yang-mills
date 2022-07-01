#ifndef YM_TRACEDEF_C
#define YM_TRACEDEF_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#ifdef OPENMP_MODE
  #include<omp.h>
#endif

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"

void real_main(char *in_file)
    {
    Gauge_Conf GC;
    Geometry geo;
    GParam param;

    char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
    int count;
    double acc, acc_local;
    FILE *datafilep, *monofilep;
    time_t time1, time2;

    // to disable nested parallelism
    #ifdef OPENMP_MODE
      omp_set_nested(0);
    #endif

    // read input file
    readinput(in_file, &param);

    // initialize random generator
    initrand(param.d_randseed);

    // open data_file
    init_data_file(&datafilep, &param);

    // open mon_file
    if(param.d_mon_meas == 1)
      {
      init_mon_file(&monofilep, &param);
      }

    // initialize geometry
    init_geometry(&geo, param.d_sizeg);

    // initialize gauge configuration
    init_gauge_conf(&GC, &geo, &param);

    // acceptance of the metropolis update
    acc=0.0;

    // montecarlo
    time(&time1);
    // count starts from 1 to avoid problems using %
    for(count=1; count < param.d_sample + 1; count++)
       {
       update_with_trace_def(&GC, &geo, &param, &acc_local);
       acc+=acc_local;

       if(count % param.d_measevery ==0 && count >= param.d_thermal)
         {
         perform_measures_localobs_with_tracedef(&GC, &geo, &param, datafilep, monofilep);
         }

       // save configuration for backup
       if(param.d_saveconf_back_every!=0)
         {
         if(count % param.d_saveconf_back_every == 0 )
           {
           // simple
           write_conf_on_file(&GC, &geo, &param);

           // backup copy
           write_conf_on_file_back(&GC, &geo, &param);
           }
         }

       // save configuration for offline analysis
       if(param.d_saveconf_analysis_every!=0)
         {
         if(count % param.d_saveconf_analysis_every == 0 )
           {
           strcpy(name, param.d_conf_file);
           sprintf(aux, "%ld", GC.update_index);
           strcat(name, aux);
           write_conf_on_file_with_name(&GC, &geo, name);
           }
         }
       }
    time(&time2);
    // montecarlo end

    acc/=(double)param.d_sample;

    // close data file
    fclose(datafilep);

    // close mon file
    if(param.d_mon_meas == 1)
      {
      fclose(monofilep);
      }

    // save configuration
    if(param.d_saveconf_back_every!=0)
      {
      write_conf_on_file(&GC, &geo, &param);
      }

    // print simulation details
    print_parameters_tracedef(&param, time1, time2, acc);

    // free gauge configuration
    free_gauge_conf(&GC, &geo);

    // free geometry
    free_geometry(&geo);
    }


void print_template_input(void)
  {
  FILE *fp;

  fp=fopen("template_input.in", "w");

  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file template_input.in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "size 4 4 4 4\n");
    fprintf(fp,"\n");
    fprintf(fp, "beta 5.705\n");
    fprintf(fp, "htracedef  1.1\n");
    fprintf(fp, "theta 1.5\n");
    fprintf(fp,"\n");
    fprintf(fp, "sample    10\n");
    fprintf(fp, "thermal   0\n");
    fprintf(fp, "overrelax 5\n");
    fprintf(fp, "measevery 1\n");
    fprintf(fp, "monomeas 0   # 1=monopoles measures are performed\n");
    fprintf(fp,"\n");
    fprintf(fp, "start                   0  # 0=ordered  1=random  2=from saved configuration\n");
    fprintf(fp, "saveconf_back_every     5  # if 0 does not save, else save backup configurations every ... updates\n");
    fprintf(fp, "saveconf_analysis_every 5  # if 0 does not save, else save configurations for analysis every ... updates\n");
    fprintf(fp, "\n");
    fprintf(fp, "epsilon_metro    0.25  #distance from the identity of the random matrix for metropolis\n");
    fprintf(fp,"\n");
    fprintf(fp, "coolsteps  3     # number of cooling steps to be used\n");
    fprintf(fp, "coolrepeat 5     # number of times 'coolsteps' are repeated\n");
    fprintf(fp,"\n");
    fprintf(fp, "#output files\n");
    fprintf(fp, "conf_file  conf.dat\n");
    fprintf(fp, "data_file  dati.dat\n");
    fprintf(fp, "mon_file   mon.dat\n");
    fprintf(fp, "log_file   log.dat\n");
    fprintf(fp, "\n");
    fprintf(fp, "randseed 0    #(0=time)\n");
    fclose(fp);
    }
  }


int main (int argc, char **argv)
    {
    char in_file[50];

    if(argc != 2)
      {
      printf("\nPackage %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
      printf("Claudio Bonati %s\n", PACKAGE_BUGREPORT);
      printf("Usage: %s input_file\n\n", argv[0]);

      printf("Compilation details:\n");
      printf("\tN_c (number of colors): %d\n", NCOLOR);
      printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
      printf("\tNum_levels (number of levels): %d\n", NLEVELS);
      printf("\n");
      printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
      printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

      #ifdef DEBUG
        printf("\n\tDEBUG mode\n");
      #endif

      #ifdef OPENMP_MODE
        printf("\n\tusing OpenMP with %d threads\n", NTHREADS);
      #endif

      #ifdef THETA_MODE
        printf("\n\tusing imaginary theta\n");
      #endif

      printf("\n");

      #ifdef __INTEL_COMPILER
        printf("\tcompiled with icc\n");
      #elif defined(__clang__)
        printf("\tcompiled with clang\n");
      #elif defined( __GNUC__ )
        printf("\tcompiled with gcc version: %d.%d.%d\n",
                __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
      #endif

      print_template_input();


      return EXIT_SUCCESS;
      }
    else
      {
      if(strlen(argv[1]) >= STD_STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in include/macro.h\n");
        }
      else
        {
        strcpy(in_file, argv[1]);
        }
      }

    real_main(in_file);

    return EXIT_SUCCESS;
    }

#endif
