#ifndef YM_POLYCORR_C
#define YM_POLYCORR_C

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
    FILE *datafilep;
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

    // initialize function_pointers
    init_function_pointers();

    // initialize geometry
    init_indexing_lexeo();
    init_geometry(&geo, &param);

    // initialize gauge configuration
    init_gauge_conf(&GC, &param);

    // initialize multilevel
    init_gauge_conf_polycorr(&GC, &param);

    // montecarlo
    time(&time1);
    // count starts from 1 to avoid problems using %
    for(count=1; count < param.d_sample + 1; count++)
       {
       update(&GC, &geo, &param);

       if(count % param.d_measevery ==0 && count >= param.d_thermal)
         {
         perform_measures_polycorr_ml(&GC, &geo, &param, datafilep);
         }

       // save configuration for backup
       if(param.d_saveconf_back_every!=0)
         {
         if(count % param.d_saveconf_back_every == 0 )
           {
           // simple
           save_on_file(&GC, &param);

           // backup copy
           save_on_file_back(&GC, &param);
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
           save_on_file_with_name(&GC, &param, name);
           }
         }
       }
    time(&time2);
    // montecarlo end

    // close data file
    fclose(datafilep);

    // save configuration
    if(param.d_saveconf_back_every!=0)
      {
      save_on_file(&GC, &param);
      }

    // print simulation details
    print_parameters_polycorr(&param, time1, time2);

    // free gauge configuration
    end_gauge_conf(&GC, &param);

    // end multilevel
    end_gauge_conf_polycorr(&GC);

    // free geometry
    free_geometry(&geo, &param);

    exit(EXIT_SUCCESS);
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
      printf("\n");
      printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
      printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

      #ifdef DEBUG
        printf("\n\tDEBUG mode\n");
      #endif

      #ifdef OPENMP_MODE
        printf("\n\tusing OpenMP with %d threads\n", NTHREADS);
      #endif

      #ifdef OPT_MULTIHIT
        printf("\tcompiled for multihit optimization\n");
      #endif

      #ifdef OPT_LEVEL1
        printf("\tcompiled for single level optimization\n");
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
