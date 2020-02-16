#ifndef GAUGE_CONF_DEF_C
#define GAUGE_CONF_DEF_C

#include"../include/macro.h"

#include<openssl/md5.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/flavour_matrix.h"
#include"../include/function_pointers.h"
#include"../include/gparam.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/tens_prod.h"
#include"../include/tens_prod_adj.h"

void init_gauge_conf(Gauge_Conf *GC, GParam const * const param)
  {
  long r, j;
  int err;

  // allocate the local lattice
  err=posix_memalign((void**) &(GC->lattice), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(GAUGE_GROUP *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**)&(GC->lattice[r]), (size_t) DOUBLE_ALIGN, (size_t )STDIM * sizeof(GAUGE_GROUP));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  #ifdef THETA_MODE
    alloc_clover_array(GC, param);
  #endif

  // initialize lattice
  if(param->d_start==0) // ordered start
    {
    GAUGE_GROUP aux1, aux2;
    one(&aux1);

    GC->update_index=0;

    for(r=0; r<(param->d_volume); r++)
       {
       for(j=0; j<STDIM; j++)
          {
          rand_matrix(&aux2);
          times_equal_real(&aux2, 0.001);
          plus_equal(&aux2, &aux1);
          unitarize(&aux2);
          equal_dag(&(GC->lattice[r][j]), &aux2);
          }
       }
    }
  if(param->d_start==1)  // random start
    {
    GAUGE_GROUP aux1;

    GC->update_index=0;

    for(r=0; r<(param->d_volume); r++)
       {
       for(j=0; j<STDIM; j++)
          {
          rand_matrix(&aux1);
          equal(&(GC->lattice[r][j]), &aux1);
          }
       }
    }

  if(param->d_start==2) // initialize from stored conf
    {
    read_gauge_conf(GC, param);
    }
  }


void read_gauge_conf(Gauge_Conf *GC, GParam const * const param)
  {
  FILE *fp;
  int i, dimension, tmp_i;
  int err, mu;
  long lex, si;
  GAUGE_GROUP matrix;
  char md5sum_new[2*MD5_DIGEST_LENGTH+1];
  char md5sum_old[2*MD5_DIGEST_LENGTH+1];

  fp=fopen(param->d_conf_file, "r"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else // read the txt header of the configuration
    {
    err=fscanf(fp, "%d", &dimension);
    if(err!=1)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(dimension != STDIM)
      {
      fprintf(stderr, "The space time dimension of the configuration (%d) does not coincide with the one of the global parameter (%d)\n",
              dimension, STDIM);
      exit(EXIT_FAILURE);
      }

    for(i=0; i<STDIM; i++)
       {
       err=fscanf(fp, "%d", &tmp_i);
       if(err!=1)
         {
         fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       if(tmp_i != param->d_size[i])
         {
         fprintf(stderr, "The size of the configuration lattice does not coincide with the one of the global parameter\n");
         exit(EXIT_FAILURE);
         }
       }

    err=fscanf(fp, "%ld %s\n", &(GC->update_index), md5sum_old);
    if(err!=2)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }

    fclose(fp);
    }

  fp=fopen(param->d_conf_file, "rb"); // open the configuration file in binary
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header
    err=0;
    while(err!='\n')
         {
         err=fgetc(fp);
         }

    for(lex=0; lex<param->d_volume; lex++)
       {
       si=lex_to_si(lex, param);
       for(mu=0; mu<STDIM; mu++)
          {
          err=read_from_binary_file_bigen(fp, &matrix);
          if(err!=0)
            {
            fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }

          equal(&(GC->lattice[si][mu]), &matrix);
          }
       }
    fclose(fp);

    // compute the new md5sum and check for consistency
    compute_md5sum_conf(md5sum_new, GC, param);
    if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
      {
      fprintf(stderr, "The computed md5sum %s does not match the stored %s for the file %s (%s, %d)\n", md5sum_new, md5sum_old, param->d_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }
  }


void free_gauge_conf(Gauge_Conf *GC, GParam const * const param)
  {
  long i;

  for(i=0; i<(param->d_volume); i++)
     {
     free(GC->lattice[i]);
     }
  free(GC->lattice);

  #ifdef THETA_MODE
    end_clover_array(GC, param);
  #endif
  }


// save a configuration in ILDG-like format
void write_conf_on_file_with_name(Gauge_Conf const * const GC,
                                  GParam const * const param,
                                  char const * const namefile)
  {
  long si, lex;
  int i, mu, err;
  char md5sum[2*MD5_DIGEST_LENGTH+1];
  FILE *fp;

  compute_md5sum_conf(md5sum, GC, param);

  fp=fopen(namefile, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%d ", STDIM);
    for(i=0; i<STDIM; i++)
       {
       fprintf(fp, "%d ", param->d_size[i]);
       }
    fprintf(fp, "%ld %s\n", GC->update_index, md5sum);
    }
  fclose(fp);

  fp=fopen(namefile, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(lex=0; lex<param->d_volume; lex++)
       {
       si=lex_to_si(lex, param);
       for(mu=0; mu<STDIM; mu++)
          {
          err=print_on_binary_file_bigen(fp, &(GC->lattice[si][mu]) );
          if(err!=0)
            {
            fprintf(stderr, "Error in writing the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    fclose(fp);
    }
  }


void write_conf_on_file(Gauge_Conf const * const GC, GParam const * const param)
  {
  write_conf_on_file_with_name(GC, param, param->d_conf_file);
  }


void write_conf_on_file_back(Gauge_Conf const * const GC, GParam const * const param)
  {
  char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
  static int counter=0;

  strcpy(name, param->d_conf_file);
  if(counter==0)
    {
    sprintf(aux, "_back0");
    }
  else
    {
    sprintf(aux, "_back1");
    }
  strcat(name, aux);

  write_conf_on_file_with_name(GC, param, name);

  counter=1-counter;
  }


// allocate GC and initialize with GC2
void init_gauge_conf_from_gauge_conf(Gauge_Conf *GC, Gauge_Conf const * const GC2, GParam const * const param) 
  {
  long r;
  int mu, err;

  // allocate the lattice
  err=posix_memalign((void**)&(GC->lattice), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(GAUGE_GROUP *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**)&(GC->lattice[r]), (size_t) DOUBLE_ALIGN, (size_t) STDIM * sizeof(GAUGE_GROUP));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  #ifdef THETA_MODE
    alloc_clover_array(GC, param);
  #endif

  // initialize GC
  for(r=0; r<(param->d_volume); r++)
     {
     for(mu=0; mu<STDIM; mu++)
        {
        equal(&(GC->lattice[r][mu]), &(GC2->lattice[r][mu]) );
        }
     }

  GC->update_index=GC2->update_index;
  }


// compute the md5sum of the configuration and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_conf(char *res, Gauge_Conf const * const GC, GParam const * const param)
  {
  MD5_CTX mdContext;
  unsigned char c[MD5_DIGEST_LENGTH];
  long si, lex;
  GAUGE_GROUP matrix;
  int mu, k;

  MD5_Init(&mdContext);
  for(lex=0; lex<param->d_volume; lex++)
     {
     si=lex_to_si(lex, param);
     for(mu=0; mu<STDIM; mu++)
        {
        equal(&matrix, &(GC->lattice[si][mu]));

        #if GGROUP == 0
          #if NCOLOR==1
            MD5_Update(&mdContext, &(matrix.comp), sizeof(double complex));
          #elif NCOLOR==2
            for(k=0; k<4; k++)
               {
               MD5_Update(&mdContext, &(matrix.comp[k]), sizeof(double));
               }
          #else
            for(k=0; k<NCOLOR*NCOLOR; k++)
               {
               MD5_Update(&mdContext, &(matrix.comp[k]), sizeof(double complex));
               }
          #endif
        #elif GGROUP == 1
          for(k=0; k<NCOLOR*NCOLOR; k++)
             {
             MD5_Update(&mdContext, &(matrix.comp[k]), sizeof(double));
             }
        #endif
        }
     }
  MD5_Final(c, &mdContext);

  for(k = 0; k < MD5_DIGEST_LENGTH; k++)
     {
     sprintf(&(res[2*k]), "%02x", c[k]);
     }
  }


// allocate the ml_polycorr arrays and related stuff
void alloc_polycorr_stuff(Gauge_Conf *GC,
                          GParam const * const param)
  {
  int i, j, err;

  err=posix_memalign((void**)&(GC->ml_polycorr), (size_t) DOUBLE_ALIGN, (size_t) NLEVELS * sizeof(TensProd **));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating ml_polycorr (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<NLEVELS; i++)
       {
       err=posix_memalign((void**)&(GC->ml_polycorr[i]), (size_t) DOUBLE_ALIGN, (size_t) (param->d_size[0] / param->d_ml_step[i]) * sizeof(TensProd *));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating ml_polycorr[%d] (%s, %d)\n", i, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       else
         {
         for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
            {
            err=posix_memalign((void**)&(GC->ml_polycorr[i][j]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(TensProd));
            if(err!=0)
              {
              fprintf(stderr, "Problems in allocating ml_polycorr[%d][%d] (%s, %d)\n", i, j, __FILE__, __LINE__);
              exit(EXIT_FAILURE);
              }
            }
         }
       }
    }

  err=posix_memalign((void**)&(GC->loc_poly), (size_t) DOUBLE_ALIGN, (size_t) (param->d_size[0]/param->d_ml_step[NLEVELS-1]) * sizeof(GAUGE_GROUP *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating loc_poly (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<param->d_size[0]/param->d_ml_step[NLEVELS-1]; i++)
       {
       err=posix_memalign((void**)&(GC->loc_poly[i]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol *sizeof(GAUGE_GROUP));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating loc_poly (%s, %d)\n", __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }
    }
  }


// free the ml_polycorr arrays and related stuff
void free_polycorr_stuff(Gauge_Conf *GC,
                         GParam const * const param)
  {
  int i, j;

  for(i=0; i<NLEVELS; i++)
     {
     for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
        {
        free(GC->ml_polycorr[i][j]);
        }
     free(GC->ml_polycorr[i]);
     }
  free(GC->ml_polycorr);

  for(i=0; i<param->d_size[0]/param->d_ml_step[NLEVELS-1]; i++)
     {
     free(GC->loc_poly[i]);
     }
  free(GC->loc_poly);
  }


// save ml_polycorr[0] arrays on file
void write_polycorr_on_file(Gauge_Conf const * const GC,
                            GParam const * const param,
                            int iteration)
  {
  long i;
  int j, err;
  char md5sum[2*MD5_DIGEST_LENGTH+1];
  FILE *fp;

  compute_md5sum_polycorr(md5sum, GC, param);

  fp=fopen(param->d_ml_file, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%ld %d %s\n", param->d_space_vol, iteration, md5sum);
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          err=print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
          if(err!=0)
            {
            fprintf(stderr, "Error in writing the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }

    fclose(fp);
    }
  }


// read ml_polycorr[0] arrays from file
void read_polycorr_from_file(Gauge_Conf const * const GC,
                             GParam const * const param,
                             int *iteration)
  {
  long i, loc_space_vol;
  int j, err;
  FILE *fp;
  char md5sum_new[2*MD5_DIGEST_LENGTH+1];
  char md5sum_old[2*MD5_DIGEST_LENGTH+1];

  fp=fopen(param->d_ml_file, "r"); // open the multilevel file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    i=fscanf(fp, "%ld %d %s\n", &loc_space_vol, iteration, md5sum_old);
    if(i!=3)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(loc_space_vol != param->d_space_vol)
      {
      fprintf(stderr, "Error: space_vol in the multilevel file %s is different from the one in the input (%s, %d)\n",
              param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "rb"); // open the multilevel file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header: loc_space_vol, iteration, hash
    i=0;
    while(i!='\n')
         {
         i=fgetc(fp);
         }

    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          err=read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
          if(err!=0)
            {
            fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }

    fclose(fp);
    }

  // compute the new md5sum and check for consistency
  compute_md5sum_polycorr(md5sum_new, GC, param);
  if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
    {
    fprintf(stderr, "The computed md5sum %s of the multilevel file does not match the stored %s\n", md5sum_new, md5sum_old);
    exit(EXIT_FAILURE);
    }
  }


// compute the md5sum of the ml_polycorr[0] arrays and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_polycorr(char *res, Gauge_Conf const * const GC, GParam const * const param)
  {
  MD5_CTX mdContext;
  unsigned char c[MD5_DIGEST_LENGTH];
  long i;
  int j;
  int n1, n2, n3, n4;

  MD5_Init(&mdContext);

  for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
     {
     for(i=0; i<(param->d_space_vol); i++)
        {
        for(n1=0; n1<NCOLOR; n1++)
           {
           for(n2=0; n2<NCOLOR; n2++)
              {
              for(n3=0; n3<NCOLOR; n3++)
                 {
                 for(n4=0; n4<NCOLOR; n4++)
                    {
                    MD5_Update(&mdContext, &((GC->ml_polycorr[0][j][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
                    }
                 }
              }
           }
        }
     }

  MD5_Final(c, &mdContext);

  for(i = 0; i < MD5_DIGEST_LENGTH; i++)
     {
     sprintf(&(res[2*i]), "%02x", c[i]);
     }
  }


// allocate the ml_polycorradj arrays
void alloc_polycorradj(Gauge_Conf *GC,
                       GParam const * const param)
  {
  int i, j, err;

  err=posix_memalign((void**)&(GC->ml_polycorradj), (size_t) DOUBLE_ALIGN, (size_t) NLEVELS * sizeof(TensProdAdj **));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating ml_polycorradj (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<NLEVELS; i++)
       {
       err=posix_memalign((void**)&(GC->ml_polycorradj[i]), (size_t) DOUBLE_ALIGN, (size_t) (param->d_size[0] / param->d_ml_step[i]) * sizeof(TensProdAdj *));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating ml_polycorradj[%d] (%s, %d)\n", i, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       else
         {
         for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
            {
            err=posix_memalign((void**)&(GC->ml_polycorradj[i][j]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(TensProdAdj));
            if(err!=0)
              {
              fprintf(stderr, "Problems in allocating ml_polycorrag[%d][%d] (%s, %d)\n", i, j, __FILE__, __LINE__);
              exit(EXIT_FAILURE);
              }
            }
         }
       }
    }

  err=posix_memalign((void**)&(GC->loc_polyadj), (size_t) DOUBLE_ALIGN, (size_t) (param->d_size[0]/param->d_ml_step[NLEVELS-1]) * sizeof(GAUGE_GROUP_ADJ *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating loc_polyadj (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<param->d_size[0]/param->d_ml_step[NLEVELS-1]; i++)
       {
       err=posix_memalign((void**)&(GC->loc_polyadj[i]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol *sizeof(GAUGE_GROUP_ADJ));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating loc_polyadj (%s, %d)\n", __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }
    }
  }


// free the ml_polycorradj arrays
void free_polycorradj(Gauge_Conf *GC,
                      GParam const * const param)
  {
  int i, j;

  for(i=0; i<NLEVELS; i++)
     {
     for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
        {
        free(GC->ml_polycorradj[i][j]);
        }
     free(GC->ml_polycorradj[i]);
     }
  free(GC->ml_polycorradj);

  for(i=0; i<param->d_size[0]/param->d_ml_step[NLEVELS-1]; i++)
     {
     free(GC->loc_polyadj[i]);
     }
  free(GC->loc_polyadj);
  }


// save ml_polycorradj[0] arrays on file
void write_polycorradj_on_file(Gauge_Conf const * const GC,
                               GParam const * const param,
                               int iteration)
  {
  long i;
  int j, err;
  char md5sum[2*MD5_DIGEST_LENGTH+1];
  FILE *fp;

  compute_md5sum_polycorradj(md5sum, GC, param);

  fp=fopen(param->d_ml_file, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%ld %d %s\n", param->d_space_vol, iteration, md5sum);
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          err=print_on_binary_file_bigen_TensProdAdj(fp, &(GC->ml_polycorradj[0][j][i]));
          if(err!=0)
            {
            fprintf(stderr, "Error in writing the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }

    fclose(fp);
    }
  }


// read ml_polycorradj[0] arrays from file
void read_polycorradj_from_file(Gauge_Conf const * const GC,
                                GParam const * const param,
                                int *iteration)
  {
  long i, loc_space_vol;
  int j, err;
  FILE *fp;
  char md5sum_new[2*MD5_DIGEST_LENGTH+1];
  char md5sum_old[2*MD5_DIGEST_LENGTH+1];

  fp=fopen(param->d_ml_file, "r"); // open the multilevel file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    i=fscanf(fp, "%ld %d %s\n", &loc_space_vol, iteration, md5sum_old);
    if(i!=3)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(loc_space_vol != param->d_space_vol)
      {
      fprintf(stderr, "Error: space_vol in the multilevel file %s is different from the one in the input (%s, %d)\n",
              param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "rb"); // open the multilevel file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header: loc_space_vol, iteration, hash
    i=0;
    while(i!='\n')
         {
         i=fgetc(fp);
         }

    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          err=read_from_binary_file_bigen_TensProdAdj(fp, &(GC->ml_polycorradj[0][j][i]));
          if(err!=0)
            {
            fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }

    fclose(fp);
    }

  // compute the new md5sum and check for consistency
  compute_md5sum_polycorradj(md5sum_new, GC, param);
  if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
    {
    fprintf(stderr, "The computed md5sum %s of the multilevel file does not match the stored %s\n", md5sum_new, md5sum_old);
    exit(EXIT_FAILURE);
    }
  }


// compute the md5sum of the ml_polycorradj[0] arrays and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_polycorradj(char *res, Gauge_Conf const * const GC, GParam const * const param)
  {
  MD5_CTX mdContext;
  unsigned char c[MD5_DIGEST_LENGTH];
  long i;
  int j;
  int n1, n2, n3, n4;

  MD5_Init(&mdContext);

  #if GGROUP == 0
    #define MAXVALUE NCOLOR*NCOLOR-1
  #elif GGROUP == 1
    #define MAXVALUE (NCOLOR*(NCOLOR-1)/2)
  #endif

  for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
     {
     for(i=0; i<(param->d_space_vol); i++)
        {
        for(n1=0; n1<MAXVALUE; n1++)
           {
           for(n2=0; n2<MAXVALUE; n2++)
              {
              for(n3=0; n3<MAXVALUE; n3++)
                 {
                 for(n4=0; n4<MAXVALUE; n4++)
                    {
                    MD5_Update(&mdContext, &((GC->ml_polycorradj[0][j][i]).comp[n1][n2][n3][n4]), sizeof(double));
                    }
                 }
              }
           }
        }
     }

  #undef MAXVALUE

  MD5_Final(c, &mdContext);

  for(i = 0; i < MD5_DIGEST_LENGTH; i++)
     {
     sprintf(&(res[2*i]), "%02x", c[i]);
     }
  }


// allocate the ml_polycorr, polyplaq arrays and related stuff
void alloc_tube_disc_stuff(Gauge_Conf *GC,
                           GParam const * const param)
  {
  int i, err;

  alloc_polycorr_stuff(GC, param);

  err=posix_memalign((void**)&(GC->ml_polyplaq), (size_t) DOUBLE_ALIGN, (size_t) NLEVELS * sizeof(TensProd *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating ml_polyplaq (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<NLEVELS; i++)
       {
       err=posix_memalign((void**)&(GC->ml_polyplaq[i]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(TensProd));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating ml_polyplaq[%d] (%s, %d)\n", i, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }
    }

  err=posix_memalign((void**)&(GC->loc_plaq), (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol *sizeof(double complex));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating loc_plaq (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// free the ml_polycorr, ml_polyplaq arrays and relates stuff
void free_tube_disc_stuff(Gauge_Conf *GC,
                          GParam const * const param)
  {
  int i;

  free_polycorr_stuff(GC, param);

  for(i=0; i<NLEVELS; i++)
     {
     free(GC->ml_polyplaq[i]);
     }
  free(GC->ml_polyplaq);

  free(GC->loc_plaq);
  }


// save ml_polycorr[0] and ml_polyplaq[0] arrays on file
void write_tube_disc_stuff_on_file(Gauge_Conf const * const GC,
                                   GParam const * const param,
                                   int iteration)
  {
  long i;
  int j, err;
  char md5sum[2*MD5_DIGEST_LENGTH+1];
  FILE *fp;

  compute_md5sum_tube_disc_stuff(md5sum, GC, param);

  fp=fopen(param->d_ml_file, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%ld %d %s\n", param->d_space_vol, iteration, md5sum);
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          err=print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
          if(err!=0)
            {
            fprintf(stderr, "Error in writing the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       err=print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaq[0][i]));
       if(err!=0)
         {
         fprintf(stderr, "Error in writing the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }

    fclose(fp);
    }
  }


// read ml_polycorr[0] and ml_polyplaq[0] arrays from file
void read_tube_disc_stuff_from_file(Gauge_Conf const * const GC,
                                    GParam const * const param,
                                    int *iteration)
  {
  long i, loc_space_vol;
  int j, err;
  FILE *fp;
  char md5sum_new[2*MD5_DIGEST_LENGTH+1];
  char md5sum_old[2*MD5_DIGEST_LENGTH+1];

  fp=fopen(param->d_ml_file, "r"); // open the multilevel file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    i=fscanf(fp, "%ld %d %s\n", &loc_space_vol, iteration, md5sum_old);
    if(i!=3)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(loc_space_vol != param->d_space_vol)
      {
      fprintf(stderr, "Error: space_vol in the multilevel file %s is different from the one in the input (%s, %d)\n",
              param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "rb"); // open the multilevel file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header: loc_space_vol, iteration, hash
    i=0;
    while(i!='\n')
         {
         i=fgetc(fp);
         }

    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          err=read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
          if(err!=0)
            {
            fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       err=read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaq[0][i]));
       if(err!=0)
         {
         fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }

    fclose(fp);
    }

  // compute the new md5sum and check for consistency
  compute_md5sum_tube_disc_stuff(md5sum_new, GC, param);
  if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
    {
    fprintf(stderr, "The computed md5sum %s of the multilevel file does not match the stored %s\n", md5sum_new, md5sum_old);
    exit(EXIT_FAILURE);
    }
  }


// compute the md5sum of the ml_polycorr[0] and ml_polyplaq[0] arrays and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_tube_disc_stuff(char *res, Gauge_Conf const * const GC, GParam const * const param)
  {
  MD5_CTX mdContext;
  unsigned char c[MD5_DIGEST_LENGTH];
  long i;
  int j;
  int n1, n2, n3, n4;

  MD5_Init(&mdContext);

  for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
     {
     for(i=0; i<(param->d_space_vol); i++)
        {
        for(n1=0; n1<NCOLOR; n1++)
           {
           for(n2=0; n2<NCOLOR; n2++)
              {
              for(n3=0; n3<NCOLOR; n3++)
                 {
                 for(n4=0; n4<NCOLOR; n4++)
                    {
                    MD5_Update(&mdContext, &((GC->ml_polycorr[0][j][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
                    }
                 }
              }
           }
        }
     }

  for(i=0; i<(param->d_space_vol); i++)
     {
     for(n1=0; n1<NCOLOR; n1++)
        {
        for(n2=0; n2<NCOLOR; n2++)
           {
           for(n3=0; n3<NCOLOR; n3++)
              {
              for(n4=0; n4<NCOLOR; n4++)
                 {
                 MD5_Update(&mdContext, &((GC->ml_polyplaq[0][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
                 }
              }
           }
        }
     }

  MD5_Final(c, &mdContext);

  for(i = 0; i < MD5_DIGEST_LENGTH; i++)
     {
     sprintf(&(res[2*i]), "%02x", c[i]);
     }
  }


// allocate the ml_polycorradj, polyplaqadj arrays and related stuff
void alloc_tubeadj_disc_stuff(Gauge_Conf *GC,
                           GParam const * const param)
  {
  int i, err;

  alloc_polycorradj(GC, param);

  err=posix_memalign((void**)&(GC->ml_polyplaqadj), (size_t) DOUBLE_ALIGN, (size_t) NLEVELS * sizeof(TensProdAdj *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating ml_polyplaqadj (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<NLEVELS; i++)
       {
       err=posix_memalign((void**)&(GC->ml_polyplaqadj[i]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(TensProdAdj));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating ml_polyplaqadj[%d] (%s, %d)\n", i, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }
    }

  err=posix_memalign((void**)&(GC->loc_plaq), (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol *sizeof(double complex));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating loc_plaq (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// free the ml_polycorradj, polyplaqadj arrays and related stuff
void free_tubeadj_disc_stuff(Gauge_Conf *GC,
                          GParam const * const param)
  {
  int i;

  free_polycorradj(GC, param);

  for(i=0; i<NLEVELS; i++)
     {
     free(GC->ml_polyplaqadj[i]);
     }
  free(GC->ml_polyplaqadj);

  free(GC->loc_plaq);
  }


// save ml_polycorradj[0], ml_polyplaqadj[0] arrays on file
void write_tubeadj_disc_stuff_on_file(Gauge_Conf const * const GC,
                                   GParam const * const param,
                                   int iteration)
  {
  long i;
  int j, err;
  char md5sum[2*MD5_DIGEST_LENGTH+1];
  FILE *fp;

  compute_md5sum_tubeadj_disc_stuff(md5sum, GC, param);

  fp=fopen(param->d_ml_file, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%ld %d %s\n", param->d_space_vol, iteration, md5sum);
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          err=print_on_binary_file_bigen_TensProdAdj(fp, &(GC->ml_polycorradj[0][j][i]));
          if(err!=0)
            {
            fprintf(stderr, "Error in writing the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       err=print_on_binary_file_bigen_TensProdAdj(fp, &(GC->ml_polyplaqadj[0][i]));
       if(err!=0)
         {
         fprintf(stderr, "Error in writing the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }

    fclose(fp);
    }
  }


// read ml_polycorradj[0]and ml_polyplaqadj[0] arrays from file
void read_tubeadj_disc_stuff_from_file(Gauge_Conf const * const GC,
                                       GParam const * const param,
                                       int *iteration)
  {
  long i, loc_space_vol;
  int j, err;
  FILE *fp;
  char md5sum_new[2*MD5_DIGEST_LENGTH+1];
  char md5sum_old[2*MD5_DIGEST_LENGTH+1];

  fp=fopen(param->d_ml_file, "r"); // open the multilevel file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    i=fscanf(fp, "%ld %d %s\n", &loc_space_vol, iteration, md5sum_old);
    if(i!=3)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(loc_space_vol != param->d_space_vol)
      {
      fprintf(stderr, "Error: space_vol in the multilevel file %s is different from the one in the input (%s, %d)\n",
              param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "rb"); // open the multilevel file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header: loc_space_vol, iteration, hash
    i=0;
    while(i!='\n')
         {
         i=fgetc(fp);
         }

    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          err=read_from_binary_file_bigen_TensProdAdj(fp, &(GC->ml_polycorradj[0][j][i]));
          if(err!=0)
            {
            fprintf(stderr, "Error in writing the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       err=read_from_binary_file_bigen_TensProdAdj(fp, &(GC->ml_polyplaqadj[0][i]));
       if(err!=0)
         {
         fprintf(stderr, "Error in writing the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }


    fclose(fp);
    }

  // compute the new md5sum and check for consistency
  compute_md5sum_tubeadj_disc_stuff(md5sum_new, GC, param);
  if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
    {
    fprintf(stderr, "The computed md5sum %s of the multilevel file does not match the stored %s\n", md5sum_new, md5sum_old);
    exit(EXIT_FAILURE);
    }
  }


// compute the md5sum of the ml_polycorradj[0] and ml_polyplaqadj[0] arrays
// and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_tubeadj_disc_stuff(char *res, Gauge_Conf const * const GC, GParam const * const param)
  {
  MD5_CTX mdContext;
  unsigned char c[MD5_DIGEST_LENGTH];
  long i;
  int j;
  int n1, n2, n3, n4;

  MD5_Init(&mdContext);

  #if GGROUP == 0
    #define MAXVALUE NCOLOR*NCOLOR-1
  #elif GGROUP == 1
    #define MAXVALUE (NCOLOR*(NCOLOR-1)/2)
  #endif

  for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
     {
     for(i=0; i<(param->d_space_vol); i++)
        {
        for(n1=0; n1<MAXVALUE; n1++)
           {
           for(n2=0; n2<MAXVALUE; n2++)
              {
              for(n3=0; n3<MAXVALUE; n3++)
                 {
                 for(n4=0; n4<MAXVALUE; n4++)
                    {
                    MD5_Update(&mdContext, &((GC->ml_polycorradj[0][j][i]).comp[n1][n2][n3][n4]), sizeof(double));
                    }
                 }
              }
           }
        }
     }

  for(i=0; i<(param->d_space_vol); i++)
     {
     for(n1=0; n1<MAXVALUE; n1++)
        {
        for(n2=0; n2<MAXVALUE; n2++)
           {
           for(n3=0; n3<MAXVALUE; n3++)
              {
              for(n4=0; n4<MAXVALUE; n4++)
                 {
                 MD5_Update(&mdContext, &((GC->ml_polyplaqadj[0][i]).comp[n1][n2][n3][n4]), sizeof(double));
                 }
              }
           }
        }
     }

  #undef MAXVALUE

  MD5_Final(c, &mdContext);

  for(i = 0; i < MD5_DIGEST_LENGTH; i++)
     {
     sprintf(&(res[2*i]), "%02x", c[i]);
     }
  }


// allocate the polycorr, polyplaq, polyplaqconn arrays and related stuff
void alloc_tube_conn_stuff(Gauge_Conf *GC,
                           GParam const * const param)
  {
  int i, err;

  alloc_tube_disc_stuff(GC, param);

  err=posix_memalign((void**)&(GC->ml_polyplaqconn), (size_t) DOUBLE_ALIGN, (size_t) NLEVELS * sizeof(TensProd *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating ml_polyplaqconn (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<NLEVELS; i++)
       {
       err=posix_memalign((void**)&(GC->ml_polyplaqconn[i]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(TensProd));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating ml_polyplaqconn[%d] (%s, %d)\n", i, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }
    }

  err=posix_memalign((void**)&(GC->loc_plaqconn), (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol *sizeof(GAUGE_GROUP));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating loc_polyplaqconn (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  }


// free the polycorr, polyplaq, polyplaqconn arrays and related stuff
void free_tube_conn_stuff(Gauge_Conf *GC,
                          GParam const * const param)
  {
  int i;

  free_tube_disc_stuff(GC, param);

  for(i=0; i<NLEVELS; i++)
     {
     free(GC->ml_polyplaqconn[i]);
     }
  free(GC->ml_polyplaqconn);

  free(GC->loc_plaqconn);
  }


// save ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] arrays on file
void write_tube_conn_stuff_on_file(Gauge_Conf const * const GC,
                                   GParam const * const param,
                                   int iteration)
  {
  long i;
  int j, err;
  char md5sum[2*MD5_DIGEST_LENGTH+1];
  FILE *fp;

  compute_md5sum_tube_conn_stuff(md5sum, GC, param);

  fp=fopen(param->d_ml_file, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%ld %d %s\n", param->d_space_vol, iteration, md5sum);
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          err=print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
          if(err!=0)
            {
            fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       err=print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaq[0][i]));
       if(err!=0)
         {
         fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       err=print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaqconn[0][i]));
       if(err!=0)
         {
         fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }

    fclose(fp);
    }
  }


// read ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] arrays from file
void read_tube_conn_stuff_from_file(Gauge_Conf const * const GC,
                                    GParam const * const param,
                                    int *iteration)
  {
  long i, loc_space_vol;
  int j, err;
  FILE *fp;
  char md5sum_new[2*MD5_DIGEST_LENGTH+1];
  char md5sum_old[2*MD5_DIGEST_LENGTH+1];

  fp=fopen(param->d_ml_file, "r"); // open the multilevel file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    i=fscanf(fp, "%ld %d %s\n", &loc_space_vol, iteration, md5sum_old);
    if(i!=3)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(loc_space_vol != param->d_space_vol)
      {
      fprintf(stderr, "Error: space_vol in the multilevel file %s is different from the one in the input (%s, %d)\n",
              param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "rb"); // open the multilevel file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header: loc_space_vol, iteration, hash
    i=0;
    while(i!='\n')
         {
         i=fgetc(fp);
         }

    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          err=read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
          if(err!=0)
            {
            fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       err=read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaq[0][i]));
       if(err!=0)
         {
         fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       err=read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaqconn[0][i]));
       if(err!=0)
         {
         fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }

    fclose(fp);
    }

  // compute the new md5sum and check for consistency
  compute_md5sum_tube_conn_stuff(md5sum_new, GC, param);
  if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
    {
    fprintf(stderr, "The computed md5sum %s of the multilevel file does not match the stored %s\n", md5sum_new, md5sum_old);
    exit(EXIT_FAILURE);
    }
  }


// compute the md5sum of the ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] arrays and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_tube_conn_stuff(char *res, Gauge_Conf const * const GC, GParam const * const param)
  {
  MD5_CTX mdContext;
  unsigned char c[MD5_DIGEST_LENGTH];
  long i;
  int j;
  int n1, n2, n3, n4;

  MD5_Init(&mdContext);

  for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
     {
     for(i=0; i<(param->d_space_vol); i++)
        {
        for(n1=0; n1<NCOLOR; n1++)
           {
           for(n2=0; n2<NCOLOR; n2++)
              {
              for(n3=0; n3<NCOLOR; n3++)
                 {
                 for(n4=0; n4<NCOLOR; n4++)
                    {
                    MD5_Update(&mdContext, &((GC->ml_polycorr[0][j][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
                    }
                 }
              }
           }
        }
     }

  for(i=0; i<(param->d_space_vol); i++)
     {
     for(n1=0; n1<NCOLOR; n1++)
        {
        for(n2=0; n2<NCOLOR; n2++)
           {
           for(n3=0; n3<NCOLOR; n3++)
              {
              for(n4=0; n4<NCOLOR; n4++)
                 {
                 MD5_Update(&mdContext, &((GC->ml_polyplaq[0][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
                 }
              }
           }
        }
     }

  for(i=0; i<(param->d_space_vol); i++)
     {
     for(n1=0; n1<NCOLOR; n1++)
        {
        for(n2=0; n2<NCOLOR; n2++)
           {
           for(n3=0; n3<NCOLOR; n3++)
              {
              for(n4=0; n4<NCOLOR; n4++)
                 {
                 MD5_Update(&mdContext, &((GC->ml_polyplaqconn[0][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
                 }
              }
           }
        }
     }

  MD5_Final(c, &mdContext);

  for(i = 0; i < MD5_DIGEST_LENGTH; i++)
     {
     sprintf(&(res[2*i]), "%02x", c[i]);
     }
  }


// allocate the clovers arrays
void alloc_clover_array(Gauge_Conf *GC,
                       GParam const * const param)
  {
  int i, err;
  long r;

  err=posix_memalign((void**)&(GC->clover_array), (size_t)DOUBLE_ALIGN, (size_t) param->d_volume *sizeof(GAUGE_GROUP **));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating clovers (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(r=0; r<param->d_volume; r++)
       {
       err=posix_memalign((void**)&(GC->clover_array[r]), (size_t)DOUBLE_ALIGN, STDIM*sizeof(GAUGE_GROUP *));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating clovers[%ld] (%s, %d)\n", r, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       for(i=0; i<STDIM; i++)
          {
          err=posix_memalign((void**)&(GC->clover_array[r][i]), (size_t)DOUBLE_ALIGN, STDIM*sizeof(GAUGE_GROUP));
          if(err!=0)
            {
            fprintf(stderr, "Problems in allocating clovers[%ld][%d] (%s, %d)\n", r, i, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    }
  }


// free the clovers arrays
void end_clover_array(Gauge_Conf *GC,
                      GParam const * const param)
  {
  int i;
  long r;

  for(r=0; r<param->d_volume; r++)
     {
     for(i=0; i<STDIM; i++)
        {
        free(GC->clover_array[r][i]);
        }
     free(GC->clover_array[r]);
     }
  free(GC->clover_array);
  }


// intialize the higgs field
void init_higgs_conf(Gauge_Conf *GC, GParam const * const param)
  {
  long r;
  int err;

  // allocate the higgs field
  err=posix_memalign((void**) &(GC->higgs), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(GAUGE_VECS));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the higgs field! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // allocate the Qh field
  err=posix_memalign((void**) &(GC->Qh), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(FMatrix));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the Qh field! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // allocate the Dh field
  err=posix_memalign((void**) &(GC->Dh), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double complex));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the Dh field! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialize the higgs field
  if(param->d_start==0) // ordered start
    {
    GAUGE_VECS aux1;

    one_vecs(&aux1);
    normalize_vecs(&aux1);

    for(r=0; r<(param->d_volume); r++)
       {
       equal_vecs(&(GC->higgs[r]), &aux1);
       }
    }

  if(param->d_start==1)  // random start
    {
    GAUGE_VECS aux1;

    for(r=0; r<(param->d_volume); r++)
       {
       rand_vecs(&aux1);
       equal_vecs(&(GC->higgs[r]), &aux1);
       }
    }

  if(param->d_start==2) // initialize from stored conf
    {
    read_higgs_conf(GC, param);
    }
  }


void read_higgs_conf(Gauge_Conf *GC, GParam const * const param)
  {
  FILE *fp;
  int i, dimension, tmp_i;
  int err;
  long lex, si;
  GAUGE_VECS vec;
  char md5sum_new[2*MD5_DIGEST_LENGTH+1];
  char md5sum_old[2*MD5_DIGEST_LENGTH+1];

  fp=fopen(param->d_higgs_conf_file, "r"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_higgs_conf_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else // read the txt header of the configuration
    {
    err=fscanf(fp, "%d", &dimension);
    if(err!=1)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_higgs_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(dimension != STDIM)
      {
      fprintf(stderr, "The space time dimension of the configuration (%d) does not coincide with the one of the global parameter (%d)\n",
              dimension, STDIM);
      exit(EXIT_FAILURE);
      }

    for(i=0; i<STDIM; i++)
       {
       err=fscanf(fp, "%d", &tmp_i);
       if(err!=1)
         {
         fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_higgs_conf_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       if(tmp_i != param->d_size[i])
         {
         fprintf(stderr, "The size of the configuration lattice does not coincide with the one of the global parameter\n");
         exit(EXIT_FAILURE);
         }
       }

    err=fscanf(fp, "%ld %s\n", &(GC->update_index), md5sum_old);
    if(err!=2)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_higgs_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }

    fclose(fp);
    }

  fp=fopen(param->d_higgs_conf_file, "rb"); // open the configuration file in binary
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header
    err=0;
    while(err!='\n')
         {
         err=fgetc(fp);
         }

    for(lex=0; lex<param->d_volume; lex++)
       {
       si=lex_to_si(lex, param);

       err=read_from_binary_file_bigen_vecs(fp, &vec);
       if(err!=0)
         {
         fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_higgs_conf_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }

       equal_vecs(&(GC->higgs[si]), &vec);
       }
    fclose(fp);

    // compute the new md5sum and check for consistency
    compute_md5sum_higgs(md5sum_new, GC, param);
    if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
      {
      fprintf(stderr, "The computed md5sum %s does not match the stored %s for the file %s (%s, %d)\n", md5sum_new, md5sum_old, param->d_higgs_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }
  }


void free_higgs_conf(Gauge_Conf *GC)
  {
  free(GC->higgs);
  free(GC->Qh);
  free(GC->Dh);
  }


// save the higgs field in ILDG-like format
void write_higgs_on_file_with_name(Gauge_Conf const * const GC,
                                   GParam const * const param,
                                   char const * const namefile)
  {
  long si, lex;
  int i, err;
  char md5sum[2*MD5_DIGEST_LENGTH+1];
  FILE *fp;

  compute_md5sum_higgs(md5sum, GC, param);

  fp=fopen(namefile, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%d ", STDIM);
    for(i=0; i<STDIM; i++)
       {
       fprintf(fp, "%d ", param->d_size[i]);
       }
    fprintf(fp, "%ld %s\n", GC->update_index, md5sum);
    }
  fclose(fp);

  fp=fopen(namefile, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(lex=0; lex<param->d_volume; lex++)
       {
       si=lex_to_si(lex, param);
       err=print_on_binary_file_bigen_vecs(fp, &(GC->higgs[si]) );
       if(err!=0)
         {
         fprintf(stderr, "Error in writing the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }
    fclose(fp);
    }
  }


void write_higgs_on_file(Gauge_Conf const * const GC, GParam const * const param)
  {
  write_higgs_on_file_with_name(GC, param, param->d_higgs_conf_file);
  }


void write_higgs_on_file_back(Gauge_Conf const * const GC, GParam const * const param)
  {
  char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
  static int counter=0;

  strcpy(name, param->d_higgs_conf_file);
  if(counter==0)
    {
    sprintf(aux, "_back0");
    }
  else
    {
    sprintf(aux, "_back1");
    }
  strcat(name, aux);

  write_higgs_on_file_with_name(GC, param, name);

  counter=1-counter;
  }


// compute the md5sum of the higgs field and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_higgs(char *res, Gauge_Conf const * const GC, GParam const * const param)
  {
  MD5_CTX mdContext;
  unsigned char c[MD5_DIGEST_LENGTH];
  long si, lex;
  GAUGE_VECS vec;
  int k;

  MD5_Init(&mdContext);
  for(lex=0; lex<param->d_volume; lex++)
     {
     si=lex_to_si(lex, param);

     equal_vecs(&vec, &(GC->higgs[si]));

     #if GGROUP==0
       #if NCOLOR==1
         for(k=0; k<NHIGGS; k++)
            {
            MD5_Update(&mdContext, &(vec.comp[k]), sizeof(double complex));
            }
       #elif NCOLOR==2
         for(k=0; k<4*NHIGGS; k++)
            {
            MD5_Update(&mdContext, &(vec.comp[k]), sizeof(double));
            }
       #else
         for(k=0; k<NHIGGS*NCOLOR; k++)
            {
            MD5_Update(&mdContext, &(vec.comp[k]), sizeof(double complex));
            }
       #endif
     #elif GGROUP==1
       for(k=0; k<NHIGGS*NCOLOR; k++)
          {
          MD5_Update(&mdContext, &(vec.comp[k]), sizeof(double));
          }
     #endif
     }
  MD5_Final(c, &mdContext);

  for(k = 0; k < MD5_DIGEST_LENGTH; k++)
     {
     sprintf(&(res[2*k]), "%02x", c[k]);
     }
  }


// alloc diagonal projection related stuff
void alloc_diag_proj_stuff(Gauge_Conf *GC,
                           GParam const * const param)
   {
   int err, dir;
   long r;

   // allocate diag_proj
   err=posix_memalign((void**) &(GC->diag_proj), (size_t)DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating diag_proj (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   for(r=0;r<param->d_volume;r++)
     {
     err=posix_memalign((void**) &(GC->diag_proj[r]), (size_t)DOUBLE_ALIGN, (size_t) STDIM * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating diag_proj (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     for(dir=0; dir<STDIM; dir++)
        {
        err=posix_memalign((void**) &(GC->diag_proj[r][dir]), (size_t)DOUBLE_ALIGN, (size_t) NCOLOR * sizeof(double));
        if(err!=0)
          {
          fprintf(stderr, "Problems in allocating diag_proj (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        }
     }

   // alloc u1_subg
   err=posix_memalign((void**) &(GC->u1_subg), (size_t)DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating u1_subg (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   for(r=0;r<param->d_volume;r++)
     {
     err=posix_memalign((void**) &(GC->u1_subg[r]), (size_t)DOUBLE_ALIGN, (size_t) STDIM * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating u1_subg (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

   // alloc uflag
   err=posix_memalign((void**) &(GC->uflag), (size_t)DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating uflag (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   for(r=0;r<param->d_volume;r++)
     {
     err=posix_memalign((void**) &(GC->uflag[r]), (size_t)DOUBLE_ALIGN, (size_t) STDIM * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating uflag (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }
   }


// free diagonal projection related stuff
void free_diag_proj_stuff(Gauge_Conf *GC,
                          GParam const * const param)
   {
   int dir;
   long r;
   
   for(r=0;r<param->d_volume;r++)
      {
      for(dir=0;dir<STDIM;dir++)
         {
         free((GC->diag_proj[r][dir]));
         }
      free((GC->diag_proj[r]));
      }
   free(GC->diag_proj);

   for(r=0;r<param->d_volume;r++)
      {
      free((GC->u1_subg[r]));
      }
   free(GC->u1_subg);

   for(r=0;r<param->d_volume;r++)
      {
      free((GC->uflag[r]));
      }
   free(GC->uflag);
   }



#endif
