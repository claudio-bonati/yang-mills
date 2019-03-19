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

  for(i=0; i<NCOLOR*NCOLOR-1; i++)
     {
     for(j=0; j<NCOLOR*NCOLOR-1; j++)
        {
        for(k=0; k<NCOLOR*NCOLOR-1; k++)
           {
           for(l=0; l<NCOLOR*NCOLOR-1; l++)
              {
              printf("%.16lf ", A->comp[i][j][k][l]);
              }
           }
        }
     }
  printf("\n");
  }


void print_on_file_TensProdAdj(FILE *fp, TensProdAdj const * const A)
  {
  int i, j, k, l, err;

  for(i=0; i<NCOLOR*NCOLOR-1; i++)
     {
     for(j=0; j<NCOLOR*NCOLOR-1; j++)
        {
        for(k=0; k<NCOLOR*NCOLOR-1; k++)
           {
           for(l=0; l<NCOLOR*NCOLOR-1; l++)
              {
              err=fprintf(fp, "%.16lf", A->comp[i][j][k][l]);
              if(err!=1)
                {
                fprintf(stderr, "Problem in writing on file a TensProdAdj (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }
              }
           }
        }
     }
  fprintf(fp, "\n");
  }


void print_on_binary_file_noswap_TensProdAdj(FILE *fp, TensProdAdj const * const A)
  {
  int i, j, k, l;
  size_t err;
  double re;

  for(i=0; i<NCOLOR*NCOLOR-1; i++)
     {
     for(j=0; j<NCOLOR*NCOLOR-1; j++)
        {
        for(k=0; k<NCOLOR*NCOLOR-1; k++)
           {
           for(l=0; l<NCOLOR*NCOLOR-1; l++)
              {
              re=A->comp[i][j][k][l];

              err=fwrite(&re, sizeof(double), 1, fp);
              if(err!=1)
                {
                fprintf(stderr, "Problem in binary writing on file a TensProdAdj (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }
              }
           }
        }
     }
  }


void print_on_binary_file_swap_TensProdAdj(FILE *fp, TensProdAdj const * const A)
  {
  int i, j, k, l;
  size_t err;
  double re;

  for(i=0; i<NCOLOR*NCOLOR-1; i++)
     {
     for(j=0; j<NCOLOR*NCOLOR-1; j++)
        {
        for(k=0; k<NCOLOR*NCOLOR-1; k++)
           {
           for(l=0; l<NCOLOR*NCOLOR-1; l++)
              {
              re=A->comp[i][j][k][l];

              SwapBytesDouble((void *)&re);

              err=fwrite(&re, sizeof(double), 1, fp);
              if(err!=1)
                {
                fprintf(stderr, "Problem in binary writing on file a TensProdAdj (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }
              }
           }
        }
     }
  }


void print_on_binary_file_bigen_TensProdAdj(FILE *fp, TensProdAdj const * const A)
  {
  if(endian()==0) // little endian machine
    {
    print_on_binary_file_swap_TensProdAdj(fp, A);
    }
  else
    {
    print_on_binary_file_noswap_TensProdAdj(fp, A);
    }
  }


void read_from_file_TensProdAdj(FILE *fp, TensProdAdj *A)
  {
  int i, j, k, l, err;
  double re;

  for(i=0; i<NCOLOR*NCOLOR-1; i++)
     {
     for(j=0; j<NCOLOR*NCOLOR-1; j++)
        {
        for(k=0; k<NCOLOR*NCOLOR-1; k++)
           {
           for(l=0; l<NCOLOR*NCOLOR-1; l++)
              {
              err=fscanf(fp, "%lg", &re);
              if(err!=1)
                {
                fprintf(stderr, "Problems reading TensProdAdj from file (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }
              A->comp[i][j][k][l]=re;
              }
           }
        }
     }
  }


void read_from_binary_file_noswap_TensProdAdj(FILE *fp, TensProdAdj *A)
  {
  int i, j, k, l;
  size_t err;
  double re;

  for(i=0; i<NCOLOR*NCOLOR-1; i++)
     {
     for(j=0; j<NCOLOR*NCOLOR-1; j++)
        {
        for(k=0; k<NCOLOR*NCOLOR-1; k++)
           {
           for(l=0; l<NCOLOR*NCOLOR-1; l++)
              {
              err=0;

              err+=fread(&re, sizeof(double), 1, fp);
              if(err!=1)
               {
               fprintf(stderr, "Problems in binary reading TensProdAdj file (%s, %d)\n", __FILE__, __LINE__);
               exit(EXIT_FAILURE);
               }

              memcpy((void *)&(A->comp[i][j][k][l]), (void*)&re, sizeof(re));
              //equivalent to A->comp[i][j][k][l]=re;
              }
           }
        }
     }
  }


void read_from_binary_file_swap_TensProdAdj(FILE *fp, TensProdAdj *A)
  {
  int i, j, k, l;
  size_t err;
  double re;

  for(i=0; i<NCOLOR*NCOLOR-1; i++)
     {
     for(j=0; j<NCOLOR*NCOLOR-1; j++)
        {
        for(k=0; k<NCOLOR*NCOLOR-1; k++)
           {
           for(l=0; l<NCOLOR*NCOLOR-1; l++)
              {
              err=0;

              err+=fread(&re, sizeof(double), 1, fp);
              if(err!=1)
               {
               fprintf(stderr, "Problems in binary reading TensProdAdj file (%s, %d)\n", __FILE__, __LINE__);
               exit(EXIT_FAILURE);
               }

              SwapBytesDouble(&re);

              memcpy((void *)&(A->comp[i][j][k][l]), (void*)&re, sizeof(re));
              //equivalent to A->comp[i][j][k][l]=re;
              }
           }
        }
     }
  }


void read_from_binary_file_bigen_TensProdAdj(FILE *fp, TensProdAdj *A)
  {
  if(endian()==0) // little endian machine
    {
    read_from_binary_file_swap_TensProdAdj(fp, A);
    }
  else
    {
    read_from_binary_file_noswap_TensProdAdj(fp, A);
    }
  }


#endif
