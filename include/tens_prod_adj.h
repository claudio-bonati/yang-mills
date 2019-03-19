#ifndef TENS_PROD_ADJ_H
#define TENS_PROD_ADJ_H

#include<complex.h>
#include<stdio.h>
#include<stdlib.h>

#include"macro.h"

// see Luscher Weisz JHEP 0109 p. 010 (2001)   (hep-lat/0108014)

typedef struct TensProdAdj {
  #if NCOLOR!=1
    double comp[NCOLOR*NCOLOR-1][NCOLOR*NCOLOR-1][NCOLOR*NCOLOR-1][NCOLOR*NCOLOR-1] __attribute__((aligned(DOUBLE_ALIGN)));
  #else // this will never be used, it is defined just to avoid warnings
     double comp[1][1][1][1] __attribute__((aligned(DOUBLE_ALIGN)));
  #endif
} TensProdAdj;


// initialize to zero
inline void zero_TensProdAdj(TensProdAdj * restrict A)
 {
 int i0, i1, i2, i3;

 for(i0=0; i0<NCOLOR*NCOLOR-1; i0++)
    {
    for(i1=0; i1<NCOLOR*NCOLOR-1; i1++)
       {
       for(i2=0; i2<NCOLOR*NCOLOR-1; i2++)
          {
          for(i3=0; i3<NCOLOR*NCOLOR-1; i3++)
             {
             A->comp[i0][i1][i2][i3]=0.0;
             }
          }
       }
    }
 }


// initialize to one
inline void one_TensProdAdj(TensProdAdj * restrict A)
 {
 int i0, i1, i2, i3;

 for(i0=0; i0<NCOLOR*NCOLOR-1; i0++)
    {
    for(i1=0; i1<NCOLOR*NCOLOR-1; i1++)
       {
       for(i2=0; i2<NCOLOR*NCOLOR-1; i2++)
          {
          for(i3=0; i3<NCOLOR*NCOLOR-1; i3++)
             {
             A->comp[i0][i1][i2][i3]=0.0;
             }
          }
       }
    }

 for(i0=0; i0<NCOLOR*NCOLOR-1; i0++)
    {
    for(i1=0; i1<NCOLOR*NCOLOR-1; i1++)
       {
       A->comp[i0][i0][i1][i1]=1.0;
       }
    }
 }


// A=B
inline void equal_TensProdAdj(TensProdAdj * restrict A, TensProdAdj const * const restrict B)
 {
 #ifdef DEBUG
 if(A==B)
   {
   fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }
 #endif

 int i0, i1, i2, i3;

 for(i0=0; i0<NCOLOR*NCOLOR-1; i0++)
    {
    for(i1=0; i1<NCOLOR*NCOLOR-1; i1++)
       {
       for(i2=0; i2<NCOLOR*NCOLOR-1; i2++)
          {
          for(i3=0; i3<NCOLOR*NCOLOR-1; i3++)
             {
             A->comp[i0][i1][i2][i3]=B->comp[i0][i1][i2][i3];
             }
          }
       }
    }
 }


// A*=r
inline void times_equal_real_TensProdAdj(TensProdAdj * restrict A, double r)
 {
 int i0, i1, i2, i3;

 for(i0=0; i0<NCOLOR*NCOLOR-1; i0++)
    {
    for(i1=0; i1<NCOLOR*NCOLOR-1; i1++)
       {
       for(i2=0; i2<NCOLOR*NCOLOR-1; i2++)
          {
          for(i3=0; i3<NCOLOR*NCOLOR-1; i3++)
             {
             A->comp[i0][i1][i2][i3]*=r;
             }
          }
       }
    }
 }


// A+=B
inline void plus_equal_TensProdAdj(TensProdAdj * restrict A, TensProdAdj const * const restrict B)
 {
 #ifdef DEBUG
 if(A==B)
   {
   fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }
 #endif

 int i0, i1, i2, i3;

 for(i0=0; i0<NCOLOR*NCOLOR-1; i0++)
    {
    for(i1=0; i1<NCOLOR*NCOLOR-1; i1++)
       {
       for(i2=0; i2<NCOLOR*NCOLOR-1; i2++)
          {
          for(i3=0; i3<NCOLOR*NCOLOR-1; i3++)
             {
             A->comp[i0][i1][i2][i3]+=B->comp[i0][i1][i2][i3];
             }
          }
       }
    }
 }


// A=B*C
inline void times_TensProdAdj(TensProdAdj * restrict A,
                               TensProdAdj const * const restrict B,
                               TensProdAdj const * const restrict C)
 {
 #ifdef DEBUG
 if(A==B || A==C || B==C)
   {
   fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }
 #endif

 int i, j, k, l;
 int a,b;

 double sum;

 for(i=0; i<NCOLOR*NCOLOR-1; i++)
    {
    for(j=0; j<NCOLOR*NCOLOR-1; j++)
       {
       for(k=0; k<NCOLOR*NCOLOR-1; k++)
          {
          for(l=0; l<NCOLOR*NCOLOR-1; l++)
             {
             sum=0.0;
             for(a=0; a<NCOLOR*NCOLOR-1; a++)
                {
                for(b=0; b<NCOLOR*NCOLOR-1; b++)
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
inline void times_equal_TensProdAdj(TensProdAdj * restrict A, TensProdAdj const * const restrict B)
 {
 #ifdef DEBUG
 if(A==B )
   {
   fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }
 #endif

 TensProdAdj tmp __attribute__((aligned(DOUBLE_ALIGN)));

 equal_TensProdAdj(&tmp, A);
 times_TensProdAdj(A, &tmp, B);
 }


inline double retr_TensProdAdj(TensProdAdj const * const restrict A)
 {
 int i0, i1;
 double tr;

 tr=0.0;
 for(i0=0; i0<NCOLOR*NCOLOR-1; i0++)
    {
    for(i1=0; i1<NCOLOR*NCOLOR-1; i1++)
       {
       tr+=A->comp[i0][i0][i1][i1];
       }
    }

 #if NCOLOR!=1    // just to avoid warnings
   tr/=((NCOLOR*NCOLOR-1)*(NCOLOR*NCOLOR-1));
 #endif

 return tr;
 }


inline double imtr_TensProdAdj(TensProdAdj const * const restrict A)
 {
 (void) A; // just to avoid warning

 return 0;
 }

void print_on_screen_TensProdAdj(TensProdAdj const * const A);
void print_on_file_TensProdAdj(FILE *fp, TensProdAdj const * const A);
void print_on_binary_file_noswap_TensProdAdj(FILE *fp, TensProdAdj const * const A);
void print_on_binary_file_swap_TensProdAdj(FILE *fp, TensProdAdj const * const A);
void print_on_binary_file_bigen_TensProdAdj(FILE *fp, TensProdAdj const * const A);
void read_from_file_TensProdAdj(FILE *fp, TensProdAdj *A);
void read_from_binary_file_noswap_TensProdAdj(FILE *fp, TensProdAdj *A);
void read_from_binary_file_swap_TensProdAdj(FILE *fp, TensProdAdj *A);
void read_from_binary_file_bigen_TensProdAdj(FILE *fp, TensProdAdj *A);


#endif // TENS_PRODAG_H



