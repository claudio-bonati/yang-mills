#ifndef TENS_PROD_C
#define TENS_PROD_C

#include<complex.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/endianness.h"
#include"../include/macro.h"
#include"../include/tens_prod.h"


void zero_TensProd(TensProd *A)
  {
  int i0, i1, i2, i3;

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        for(i2=0; i2<NCOLOR; i2++)
           {
           for(i3=0; i3<NCOLOR; i3++)
              {
              A->comp[i0][i1][i2][i3]=0.0+0.0*I;
              }
           }
        }
     }
  }


void one_TensProd(TensProd *A)
  {
  int i0, i1, i2, i3;

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        for(i2=0; i2<NCOLOR; i2++)
           {
           for(i3=0; i3<NCOLOR; i3++)
              {
              A->comp[i0][i1][i2][i3]=0.0+0.0*I;
              }
           }
        }
     }

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        A->comp[i0][i0][i1][i1]=1.0+0.0*I;
        }
     }
  }


// A=B
void equal_TensProd(TensProd * restrict A, TensProd const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
  A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
  B  = __builtin_assume_aligned(B, DOUBLE_ALIGN);
  #else
    #ifdef __INTEL_COMPILER
     __assume_aligned(A, DOUBLE_ALIGN);
     __assume_aligned(B, DOUBLE_ALIGN);
    #endif
  #endif

  int i0, i1, i2, i3;

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        for(i2=0; i2<NCOLOR; i2++)
           {
           for(i3=0; i3<NCOLOR; i3++)
              {
              A->comp[i0][i1][i2][i3]=B->comp[i0][i1][i2][i3];
              }
           }
        }
     }
  }


// A*=r real
void times_equal_real_TensProd(TensProd * restrict A, double r)
  {
  #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
  A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
  #else
    #ifdef __INTEL_COMPILER
     __assume_aligned(A, DOUBLE_ALIGN);
    #endif
  #endif

  int i0, i1, i2, i3;

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        for(i2=0; i2<NCOLOR; i2++)
           {
           for(i3=0; i3<NCOLOR; i3++)
              {
              A->comp[i0][i1][i2][i3]*=r;
              }
           }
        }
     }
  }


// A*=r complex
void times_equal_comples_TensProd(TensProd * restrict A, double complex r)
  {
  #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
  A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
  #else
    #ifdef __INTEL_COMPILER
     __assume_aligned(A, DOUBLE_ALIGN);
    #endif
  #endif

  int i0, i1, i2, i3;

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        for(i2=0; i2<NCOLOR; i2++)
           {
           for(i3=0; i3<NCOLOR; i3++)
              {
              A->comp[i0][i1][i2][i3]*=r;
              }
           }
        }
     }
  }


// A+=B
void plus_equal_TensProd(TensProd * restrict A, TensProd const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
  A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
  B  = __builtin_assume_aligned(B, DOUBLE_ALIGN);
  #else
    #ifdef __INTEL_COMPILER
     __assume_aligned(A, DOUBLE_ALIGN);
     __assume_aligned(B, DOUBLE_ALIGN);
    #endif
  #endif

  int i0, i1, i2, i3;

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        for(i2=0; i2<NCOLOR; i2++)
           {
           for(i3=0; i3<NCOLOR; i3++)
              {
              A->comp[i0][i1][i2][i3]+=B->comp[i0][i1][i2][i3];
              }
           }
        }
     }
  }


// A=B*C
void times_TensProd(TensProd * restrict A,
                    TensProd const * restrict  B,
                    TensProd const * restrict  C)
  {
  #ifdef DEBUG
  if(A==B || A==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
  A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
  B  = __builtin_assume_aligned(B, DOUBLE_ALIGN);
  C  = __builtin_assume_aligned(C, DOUBLE_ALIGN);
  #else
    #ifdef __INTEL_COMPILER
     __assume_aligned(A, DOUBLE_ALIGN);
     __assume_aligned(B, DOUBLE_ALIGN);
     __assume_aligned(C, DOUBLE_ALIGN);
    #endif
  #endif

  int i, j, k, l;
  int a,b;

  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              sum=0.0+0.0*I;
              for(a=0; a<NCOLOR; a++)
                 {
                 for(b=0; b<NCOLOR; b++)
                    {
                    sum+=B->comp[i][a][k][b] * C->comp[a][j][b][l];
                    }
                 }
              A->comp[i][j][k][l]=sum;
              }
           }
        }
     }
  }


// A*=B
void times_equal_TensProd(TensProd * restrict A, TensProd const * restrict B)
  {
  #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
  A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
  B  = __builtin_assume_aligned(B, DOUBLE_ALIGN);
  #else
    #ifdef __INTEL_COMPILER
     __assume_aligned(A, DOUBLE_ALIGN);
     __assume_aligned(B, DOUBLE_ALIGN);
    #endif
  #endif

  TensProd tmp __attribute__((aligned(DOUBLE_ALIGN)));

  equal_TensProd(&tmp, A);
  times_TensProd(A, &tmp, B);
  }


double retr_TensProd(TensProd const * restrict A)
  {
  #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
  A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
  #else
    #ifdef __INTEL_COMPILER
     __assume_aligned(A, DOUBLE_ALIGN);
    #endif
  #endif

  int i0, i1;
  double complex tr;
  double ris;

  tr=0.0;
  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        tr+=A->comp[i0][i0][i1][i1];
        }
     }
  ris=creal(tr);
  ris/=(NCOLOR*NCOLOR);

  return ris;
  }


double imtr_TensProd(TensProd const * restrict A)
  {
  #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
  A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
  #else
    #ifdef __INTEL_COMPILER
     __assume_aligned(A, DOUBLE_ALIGN);
    #endif
  #endif

  int i0, i1;
  double complex tr;
  double ris;

  tr=0.0;
  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        tr+=A->comp[i0][i0][i1][i1];
        }
     }
  ris=cimag(tr);
  ris/=(NCOLOR*NCOLOR);

  return ris;
  }


void print_on_screen_TensProd(TensProd const * const A)
  {
  int i, j, k, l;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              printf("%.16lf %.16lf ", creal(A->comp[i][j][k][l]), cimag(A->comp[i][j][k][l]) );
              }
           }
        }
     }
  printf("\n");
  }


void print_on_file_TensProd(FILE *fp, TensProd const * const A)
  {
  int i, j, k, l, err;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              err=fprintf(fp, "%.16lf %.16lf", creal(A->comp[i][j][k][l]), cimag(A->comp[i][j][k][l]) );
              if(err!=2)
                {
                fprintf(stderr, "Problem in writing on file a TensProd (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }
              }
           }
        }
     }
  fprintf(fp, "\n");
  }


void print_on_binary_file_noswap_TensProd(FILE *fp, TensProd const * const A)
  {
  int i, j, k, l;
  size_t err;
  double re, im;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              re=creal(A->comp[i][j][k][l]);
              im=cimag(A->comp[i][j][k][l]);

              err=fwrite(&re, sizeof(double), 1, fp);
              if(err!=1)
                {
                fprintf(stderr, "Problem in binary writing on file a TensProd (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }
              err=fwrite(&im, sizeof(double), 1, fp);
              if(err!=1)
                {
                fprintf(stderr, "Problem in binary writing on file a TensProd (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }
              }
           }
        }
     }
  }


void print_on_binary_file_swap_TensProd(FILE *fp, TensProd const * const A)
  {
  int i, j, k, l;
  size_t err;
  double re, im;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              re=creal(A->comp[i][j][k][l]);
              im=cimag(A->comp[i][j][k][l]);

              SwapBytesDouble((void *)&re);
              SwapBytesDouble((void *)&im);

              err=fwrite(&re, sizeof(double), 1, fp);
              if(err!=1)
                {
                fprintf(stderr, "Problem in binary writing on file a TensProd (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }
              err=fwrite(&im, sizeof(double), 1, fp);
              if(err!=1)
                {
                fprintf(stderr, "Problem in binary writing on file a TensProd (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }
              }
           }
        }
     }
  }


void print_on_binary_file_bigen_TensProd(FILE *fp, TensProd const * const A)
  {
  if(endian()==0) // little endian machine
    {
    print_on_binary_file_swap_TensProd(fp, A);
    }
  else
    {
    print_on_binary_file_noswap_TensProd(fp, A);
    }
  }


void read_from_file_TensProd(FILE *fp, TensProd *A)
  {
  int i, j, k, l, err;
  double re, im;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              err=fscanf(fp, "%lg %lg", &re, &im);
              if(err!=2)
                {
                fprintf(stderr, "Problems reading TensProd from file (%s, %d)\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
                }
              A->comp[i][j][k][l]=re+I*im;
              }
           }
        }
     }
  }


void read_from_binary_file_noswap_TensProd(FILE *fp, TensProd *A)
  {
  int i, j, k, l;
  size_t err;
  double re, im;
  double aux[2];

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              err=0;

              err+=fread(&re, sizeof(double), 1, fp);
              err+=fread(&im, sizeof(double), 1, fp);
              if(err!=2)
               {
               fprintf(stderr, "Problems in binary reading TensProd file (%s, %d)\n", __FILE__, __LINE__);
               exit(EXIT_FAILURE);
               }

              aux[0]=re;
              aux[1]=im;

              memcpy((void *)&(A->comp[i][j][k][l]), (void*)aux, sizeof(aux));
              //equivalent to A->comp[i][j][k][l]=re+im*I;
              }
           }
        }
     }
  }


void read_from_binary_file_swap_TensProd(FILE *fp, TensProd *A)
  {
  int i, j, k, l;
  size_t err;
  double re, im;
  double aux[2];

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              err=0;

              err+=fread(&re, sizeof(double), 1, fp);
              err+=fread(&im, sizeof(double), 1, fp);
              if(err!=2)
               {
               fprintf(stderr, "Problems in binary reading TensProd file (%s, %d)\n", __FILE__, __LINE__);
               exit(EXIT_FAILURE);
               }

              SwapBytesDouble(&re);
              SwapBytesDouble(&im);
              aux[0]=re;
              aux[1]=im;

              memcpy((void *)&(A->comp[i][j][k][l]), (void*)aux, sizeof(aux));
              //equivalent to A->comp[i][j][k][l]=re+im*I;
              }
           }
        }
     }
  }


void read_from_binary_file_bigen_TensProd(FILE *fp, TensProd *A)
  {
  if(endian()==0) // little endian machine
    {
    read_from_binary_file_swap_TensProd(fp, A);
    }
  else
    {
    read_from_binary_file_noswap_TensProd(fp, A);
    }
  }

#endif
