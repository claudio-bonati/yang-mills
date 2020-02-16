#ifndef TENS_PROD_ADJ_C
#define TENS_PROD_ADJ_C

#include<complex.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/aligncheck.h"
#include"../include/endianness.h"
#include"../include/macro.h"
#include"../include/tens_prod_adj.h"

#if GGROUP == 0
  #define MAXINDEX NCOLOR*NCOLOR -1
#elif GGROUP == 1
  #define MAXINDEX NCOLOR*(NCOLOR-1)/2
#endif

// initialize to zero
void zero_TensProdAdj(TensProdAdj * restrict A);


// initialize to one
void one_TensProdAdj(TensProdAdj * restrict A);


// A=B
void equal_TensProdAdj(TensProdAdj * restrict A, TensProdAdj const * const restrict B);


// A*=r real
void times_equal_real_TensProdAdj(TensProdAdj * restrict A, double r);

// A+=B
void plus_equal_TensProdAdj(TensProdAdj * restrict A, TensProdAdj const * const restrict B);


// A=B*C
void times_TensProdAdj(TensProdAdj * restrict A,
                        TensProdAdj const * const restrict B,
                        TensProdAdj const * const restrict C);


// A*=B
void times_equal_TensProdAdj(TensProdAdj * restrict A, TensProdAdj const * const restrict B);


double retr_TensProdAdj(TensProdAdj const * const restrict A);
double imtr_TensProdAdj(TensProdAdj const * const restrict A);


void print_on_screen_TensProdAdj(TensProdAdj const * const A)
  {
  int i, j, k, l;

  for(i=0; i<MAXINDEX; i++)
     {
     for(j=0; j<MAXINDEX; j++)
        {
        for(k=0; k<MAXINDEX; k++)
           {
           for(l=0; l<MAXINDEX; l++)
              {
              printf("%.16lf ", A->comp[i][j][k][l]);
              }
           }
        }
     }
  printf("\n");
  }


int print_on_file_TensProdAdj(FILE *fp, TensProdAdj const * const A)
  {
  int i, j, k, l, err;

  for(i=0; i<MAXINDEX; i++)
     {
     for(j=0; j<MAXINDEX; j++)
        {
        for(k=0; k<MAXINDEX; k++)
           {
           for(l=0; l<MAXINDEX; l++)
              {
              err=fprintf(fp, "%.16lf", A->comp[i][j][k][l]);
              if(err!=1)
                {
                fprintf(stderr, "Problem in writing on file a TensProdAdj (%s, %d)\n", __FILE__, __LINE__);
                return 1;
                }
              }
           }
        }
     }
  fprintf(fp, "\n");
  return 0;
  }


int print_on_binary_file_noswap_TensProdAdj(FILE *fp, TensProdAdj const * const A)
  {
  int i, j, k, l;
  size_t err;
  double re;

  for(i=0; i<MAXINDEX; i++)
     {
     for(j=0; j<MAXINDEX; j++)
        {
        for(k=0; k<MAXINDEX; k++)
           {
           for(l=0; l<MAXINDEX; l++)
              {
              re=A->comp[i][j][k][l];

              err=fwrite(&re, sizeof(double), 1, fp);
              if(err!=1)
                {
                fprintf(stderr, "Problem in binary writing on file a TensProdAdj (%s, %d)\n", __FILE__, __LINE__);
                return 1;
                }
              }
           }
        }
     }
  return 0;
  }


int print_on_binary_file_swap_TensProdAdj(FILE *fp, TensProdAdj const * const A)
  {
  int i, j, k, l;
  size_t err;
  double re;

  for(i=0; i<MAXINDEX; i++)
     {
     for(j=0; j<MAXINDEX; j++)
        {
        for(k=0; k<MAXINDEX; k++)
           {
           for(l=0; l<MAXINDEX; l++)
              {
              re=A->comp[i][j][k][l];

              SwapBytesDouble((void *)&re);

              err=fwrite(&re, sizeof(double), 1, fp);
              if(err!=1)
                {
                fprintf(stderr, "Problem in binary writing on file a TensProdAdj (%s, %d)\n", __FILE__, __LINE__);
                return 1;
                }
              }
           }
        }
     }

  return 0;
  }


int print_on_binary_file_bigen_TensProdAdj(FILE *fp, TensProdAdj const * const A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=print_on_binary_file_swap_TensProdAdj(fp, A);
    }
  else
    {
    err=print_on_binary_file_noswap_TensProdAdj(fp, A);
    }

  return err;
  }


int read_from_file_TensProdAdj(FILE *fp, TensProdAdj *A)
  {
  int i, j, k, l, err;
  double re;

  for(i=0; i<MAXINDEX; i++)
     {
     for(j=0; j<MAXINDEX; j++)
        {
        for(k=0; k<MAXINDEX; k++)
           {
           for(l=0; l<MAXINDEX; l++)
              {
              err=fscanf(fp, "%lg", &re);
              if(err!=1)
                {
                fprintf(stderr, "Problems reading TensProdAdj from file (%s, %d)\n", __FILE__, __LINE__);
                return 1;
                }
              A->comp[i][j][k][l]=re;
              }
           }
        }
     }

  return 0;
  }


int read_from_binary_file_noswap_TensProdAdj(FILE *fp, TensProdAdj *A)
  {
  int i, j, k, l;
  size_t err;
  double re;

  for(i=0; i<MAXINDEX; i++)
     {
     for(j=0; j<MAXINDEX; j++)
        {
        for(k=0; k<MAXINDEX; k++)
           {
           for(l=0; l<MAXINDEX; l++)
              {
              err=0;

              err+=fread(&re, sizeof(double), 1, fp);
              if(err!=1)
               {
               fprintf(stderr, "Problems in binary reading TensProdAdj file (%s, %d)\n", __FILE__, __LINE__);
               return 1;
               }

              memcpy((void *)&(A->comp[i][j][k][l]), (void*)&re, sizeof(re));
              //equivalent to A->comp[i][j][k][l]=re;
              }
           }
        }
     }
  return 0;
  }


int read_from_binary_file_swap_TensProdAdj(FILE *fp, TensProdAdj *A)
  {
  int i, j, k, l;
  size_t err;
  double re;

  for(i=0; i<MAXINDEX; i++)
     {
     for(j=0; j<MAXINDEX; j++)
        {
        for(k=0; k<MAXINDEX; k++)
           {
           for(l=0; l<MAXINDEX; l++)
              {
              err=0;

              err+=fread(&re, sizeof(double), 1, fp);
              if(err!=1)
               {
               fprintf(stderr, "Problems in binary reading TensProdAdj file (%s, %d)\n", __FILE__, __LINE__);
               return 1;
               }

              SwapBytesDouble(&re);

              memcpy((void *)&(A->comp[i][j][k][l]), (void*)&re, sizeof(re));
              //equivalent to A->comp[i][j][k][l]=re;
              }
           }
        }
     }
  return 0;
  }


int read_from_binary_file_bigen_TensProdAdj(FILE *fp, TensProdAdj *A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=read_from_binary_file_swap_TensProdAdj(fp, A);
    }
  else
    {
    err=read_from_binary_file_noswap_TensProdAdj(fp, A);
    }

  return err;
  }

#undef MAXINDEX

#endif
