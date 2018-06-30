#ifndef YM_LOCAL_C
#define YM_LOCAL_C

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
    Gauge_Conf GC, help1, help2;
    Geometry geo;
    GParam param;

    long count1, count2;
    const long max_count=1000;
    double gftime, energy;

    FILE *datafilep;
    time_t time1, time2;

    // to disable nested parallelism
    #ifdef OPENMP_MODE
      omp_set_nested(0);
    #endif

    // read input file
    readinput(in_file, &param);

    // this code has to start from saved conf.
    param.d_start=2;

    // initialize random generator
    initrand(param.d_randseed);

    // open data_file
    init_data_file(&datafilep, &param);

    // initialize function_pointers
    init_function_pointers();

    // initialize geometry
    init_indexing_lexeo(&geo);
    init_geometry(&geo, &param);

    // initialize gauge configurations
    init_gauge_conf(&GC, &param);
    init_gauge_conf_from_gauge_conf(&help1, &GC, &param);
    init_gauge_conf_from_gauge_conf(&help2, &GC, &param);

    time(&time1);
    gftime=0.0;
    count1=0;
    while(count1<max_count)
         {
         for(count2=0; count2<20; count2++)
            {
            gradflow_RKstep(&GC, &help1, &help2, &geo, &param, param.d_gfstep);
            gftime+=param.d_gfstep;
            }

         clover_disc_energy(&GC, &geo, &param, &energy);
         printf("%lf  %lf\n", energy, energy*gftime*gftime);
         fflush(stdout);
         count1++;
         if(energy*gftime*gftime>0.35)
           {
           count1+=(max_count+10);
           }
         }
    if(count1==max_count)
      {
      fprintf(stderr, "max_count reached in (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    time(&time2);

    // close data file
    fclose(datafilep);

    // print simulation details
    print_parameters_t0(&param, time1, time2);

    // free gauge configurations
    end_gauge_conf(&GC, &param);
    end_gauge_conf(&help1, &param);
    end_gauge_conf(&help2, &param);

    // free geometry
    free_geometry(&geo, &param);

    exit(EXIT_SUCCESS);
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
    fprintf(fp, "#for gradient flow evolution\n");
    fprintf(fp, "gfstep   0.01    # integration step for gradient flow\n");
    fprintf(fp, "\n");
    fprintf(fp, "#output files\n");
    fprintf(fp, "conf_file  conf.dat\n");
    fprintf(fp, "data_file  dati.dat\n");
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

