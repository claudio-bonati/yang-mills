#ifndef YM_POLYCORR_HIGGS_LONG_C
#define YM_POLYCORR_HIGGS_LONG_C

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

    int count;
    double acc, acc_local;
    FILE *datafilep;
    time_t time1, time2;

    // to disable nested parallelism
    #ifdef OPENMP_MODE
      omp_set_nested(0);
    #endif

    // read input file
    readinput(in_file, &param);

    int tmp=param.d_sizeg[1];
    for(count=2; count<STDIM; count++)
       {
       if(tmp!= param.d_sizeg[count])
         {
         fprintf(stderr, "When using yang_mills_polycorr_higgs_long all the spatial sizes have to be of equal length.\n");
         exit(EXIT_FAILURE);
         }
       }

    // initialize random generator
    initrand(param.d_randseed);

    // open data_file
    init_data_file(&datafilep, &param);

    // initialize geometry
    init_geometry(&geo, param.d_sizeg);

    // initialize gauge configuration
    init_gauge_conf(&GC, &geo, &param);
    init_higgs_conf(&GC, &geo, &param);

    // initialize ml_polycorr arrays
    alloc_polycorr_stuff(&GC, &geo, &param);

    // acceptance of the metropolis update
    acc=0.0;

    // montecarlo starts
    time(&time1);
    if(param.d_start != 2) // NEW SIMULATION
      {
      for(count=0; count<param.d_measevery; count++)
         {
         update_with_higgs(&GC, &geo, &param, &acc_local);
         }

      if(count>param.d_thermal)
        {
        acc+=acc_local;
        }

      // save configuration
      write_conf_on_file(&GC, &geo, &param);
      write_higgs_on_file(&GC, &geo, &param);

      // backup copy
      write_conf_on_file_back(&GC, &geo, &param);
      write_higgs_on_file_back(&GC, &geo, &param);

      // save ml polycorr arrays
      write_polycorr_on_file(&GC, &geo, &param, 0);
      }
    else // CONTINUATION OF PREVIOUS SIMULATION
      {
      int count, iteration;

      // read multilevel stuff
      read_polycorr_from_file(&GC, &geo, &param, &iteration);

      if(iteration<0) // update the conf, no multilevel
        {
        for(count=0; count<param.d_measevery; count++)
           {
           update_with_higgs(&GC, &geo, &param, &acc_local);
           }

        if(count>param.d_thermal)
          {
          acc+=acc_local;
          }

        // save configuration
        write_conf_on_file(&GC, &geo, &param);
        write_higgs_on_file(&GC, &geo, &param);

        // backup copy
        write_conf_on_file_back(&GC, &geo, &param);
        write_higgs_on_file(&GC, &geo, &param);

        // save multilevel stuff
        write_polycorr_on_file(&GC, &geo, &param, 0);
        }
      else // iteration >=0, perform multilevel
        {
        multilevel_polycorr_long_with_higgs(&GC,
                                            &geo,
                                            &param,
                                            param.d_ml_step[0],
                                            iteration);
        iteration+=1;
        if(iteration==param.d_ml_level0_repeat)
          {
          // print the measure
          perform_measures_polycorr_long(&GC, &geo, &param, datafilep);

          iteration=-1; // next time the conf will be updated, no multilevel
          }

        // save multilevel stuff
        write_polycorr_on_file(&GC, &geo, &param, iteration);
        }
      }
    time(&time2);
    // montecarlo end

    acc/=(double)(param.d_sample-param.d_thermal);

    // close data file
    fclose(datafilep);

    // save configuration
    if(param.d_saveconf_back_every!=0)
      {
      write_conf_on_file(&GC, &geo, &param);
      write_higgs_on_file(&GC, &geo, &param);
      }

    // print simulation details
    print_parameters_polycorr_higgs_long(&param, time1, time2, acc);

    // free gauge configuration
    free_gauge_conf(&GC, &geo);
    free_higgs_conf(&GC);

    // free ml_polycorr
    free_polycorr_stuff(&GC, &geo, &param);

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
    fprintf(fp, "higgs_beta 1.5\n");
    fprintf(fp, "theta 1.5\n");
    fprintf(fp,"\n");
    fprintf(fp, "overrelax 5\n");
    fprintf(fp, "measevery 1\n");
    fprintf(fp,"\n");
    fprintf(fp, "start                   0  # 0=ordered  1=random  2=from saved configuration\n");
    fprintf(fp,"\n");
    fprintf(fp, "#for multilevel\n");
    fprintf(fp, "multihit         10  # number of multihit step\n");
    fprintf(fp, "ml_step          2   # timeslices for multilevel (from largest to smallest)\n");
    fprintf(fp, "ml_upd           10  # number of updates for various levels\n");
    fprintf(fp, "ml_level0_repeat 1   # number of times level0 is repeated in long sim.\n");
    fprintf(fp, "dist_poly        2   # distance between the polyakov loop\n");
    fprintf(fp,"\n");
    fprintf(fp, "epsilon_metro    0.25      #distance from the identity of the random matrix for metropolis\n");
    fprintf(fp,"\n");
    fprintf(fp, "#output files\n");
    fprintf(fp, "conf_file  conf.dat\n");
    fprintf(fp, "higgs_conf_file  higgs_conf.dat\n");
    fprintf(fp, "data_file  dati.dat\n");
    fprintf(fp, "log_file   log.dat\n");
    fprintf(fp, "ml_file    ml.dat\n");
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
      printf("\tN_higgs (number of higgs flavours): %d\n", NHIGGS);
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

      #ifdef OPT_MULTILEVEL
        printf("\tcompiled for multilevel optimization\n");
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

