#ifndef GPARAM_C
#define GPARAM_C

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/endianness.h"
#include"../include/gparam.h"
#include"../include/macro.h"

void readinput(char *in_file, GParam *param)
    {
    FILE *input;
    char str[STD_STRING_LENGTH], temp_str[STD_STRING_LENGTH];
    char c;
    double temp_d;
    int temp_i, i;
    int err, end=1;
    unsigned int temp_ui;
    
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
           err=fscanf(input, "%s", str);
           if(err!=1)
             {
             fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
             printf("err=%d", err);
             exit(EXIT_FAILURE);
             }

           // if there is match, copy the value and the next carachter,
           // to check if it is \n or if there are some comments
           // see at the end of the "else if"

           if(strncmp(str, "size", 4)==0)
             {
             for(i=0; i<STDIM; i++)
                {
                err=fscanf(input, "%d%c", &temp_i, &c);
                if(err!=2)
                  {
                  fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                  exit(EXIT_FAILURE);
                  }
                param->d_size[i]=temp_i;
                }
             }

           else if(strncmp(str, "beta", 4)==0)
                  { 
                  err=fscanf(input, "%lf%c", &temp_d, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_beta=temp_d;
                  }

           else if(strncmp(str, "sample", 6)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_sample=temp_i;
                  }
           else if(strncmp(str, "thermal", 7)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_thermal=temp_i;
                  }
           else if(strncmp(str, "overrelax", 9)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_overrelax=temp_i;
                  }
           else if(strncmp(str, "measevery", 9)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_measevery=temp_i;
                  }

           else if(strncmp(str, "start", 5)==0)
                  { 
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_start=temp_i;
                  }
           else if(strncmp(str, "saveconf_back_every", 19)==0)
                  {
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_saveconf_back_every=temp_i;
                  }
           else if(strncmp(str, "saveconf_analysis_every", 23)==0)
                  {
                  err=fscanf(input, "%d%c", &temp_i, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_saveconf_analysis_every=temp_i;
                  }

           else if(strncmp(str, "conf_file", 9)==0)
                  { 
                  err=fscanf(input, "%s%c", temp_str, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_conf_file, temp_str);
                  }
           else if(strncmp(str, "data_file", 9)==0)
                  { 
                  err=fscanf(input, "%s%c", temp_str, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_data_file, temp_str);
                  }
           else if(strncmp(str, "log_file", 8)==0)
                  { 
                  err=fscanf(input, "%s%c", temp_str, &c);
                  if(err!=2)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_log_file, temp_str);
                  }

           else if(strncmp(str, "randseed", 8)==0)
                  { 
                  err=fscanf(input, "%u%c", &temp_ui, &c);
                  if(err!=2)
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
   
           // discard eventual comments
           if(c!='\n')
             { 
             do
               {
               temp_i=getc(input);
               }
             while(temp_i!='\n');
             }

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
      #ifdef OPENMP_MODE
      for(i=0; i<STDIM; i++)
         {
         temp_i = param->d_size[i] % 2;
         if(temp_i!=0)
           {
           fprintf(stderr, "ERROR: size[%d] is not even.\n", i);
           fprintf(stderr, "When using OpenMP all the sides of the lattice have to be even! (%s, %d)\n", __FILE__, __LINE__);
           exit(EXIT_FAILURE);
           }
         }
      #endif

      init_derived_constants(param);
      }
    }


void init_derived_constants(GParam *param)
  {
  int i;

  // derived constants
  param->d_volume=1;
  for(i=0; i<STDIM; i++)
     {
     (param->d_volume)*=(param->d_size[i]);
     }
  param->d_space_vol=1;

  // direction 0 is time
  for(i=1; i<STDIM; i++)
     {
     (param->d_space_vol)*=(param->d_size[i]);
     }
  param->d_inv_vol=1.0/((double) param->d_volume);
  param->d_inv_space_vol=1.0/((double) param->d_space_vol);
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
      fprintf(*dataf, "%s ", QUOTEME(STDIM));
      for(i=0; i<STDIM; i++)
         {
         fprintf(*dataf, "%d ", param->d_size[i]);
         }
      fprintf(*dataf, "\n");
      }
    }
  else
    {
    *dataf=fopen(param->d_data_file, "w");
    fprintf(*dataf, "%s ", QUOTEME(STDIM));
    for(i=0; i<STDIM; i++)
       {
       fprintf(*dataf, "%d ", param->d_size[i]);
       }
    fprintf(*dataf, "\n");
    }
  fflush(*dataf);
  }


// print simulation parameters
void print_parameters(GParam const * const param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+--------------------+\n");
    fprintf(fp, "| Simulation details |\n");
    fprintf(fp, "+--------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %s\n\n", QUOTEME(STDIM));

    fprintf(fp, "lattice: %d", param->d_size[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_size[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
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

