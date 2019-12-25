#ifndef FLAVOUR_MATRIX_H
#define FLAVOUR_MATRIX_H

#include<complex.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"macro.h"

typedef struct FMatrix {
   double complex comp[NHIGGS*NHIGGS] __attribute__((aligned(DOUBLE_ALIGN)));
} FMatrix;
//
//  the element [i][j] can be obtained by matrix.comp[mf(i,j)] with mf(i,j) defined in macro.h
//



// A=0
inline void zero_FMatrix(FMatrix * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NHIGGS*NHIGGS; i++)
     {
     A->comp[i]=0.0+0.0*I;
     }
  }


// A=B
inline void equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NHIGGS*NHIGGS; i++)
     {
     A->comp[i]=B->comp[i];
     }
  }


// A+=B
inline void plus_equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NHIGGS*NHIGGS; i++)
     {
     A->comp[i]+=B->comp[i];
     }
  }


// A-=B
inline void minus_equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NHIGGS*NHIGGS; i++)
     {
     A->comp[i]-=B->comp[i];
     }
  }


// A-=B^{dag}
inline void minus_equal_dag_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<NHIGGS; j++)
        {
        A->comp[mf(i,j)]-=conj(B->comp[mf(j,i)]);
        }
     }
  }


// A*=r
inline void times_equal_real_FMatrix(FMatrix * restrict A, double r)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NHIGGS*NHIGGS; i++)
     {
     A->comp[i]*=r;
     }
  }


// A*=r
inline void times_equal_complex_FMatrix(FMatrix * restrict A, double complex r)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NHIGGS*NHIGGS; i++)
     {
     A->comp[i]*=r;
     }
  }


// A*=B
inline void times_equal_FMatrix(FMatrix * restrict A, FMatrix const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
   {
   fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex aux[NHIGGS] __attribute__((aligned(DOUBLE_ALIGN)));
  double complex sum;

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<NHIGGS; j++)
        {
        aux[j]=A->comp[mf(i,j)];
        }

     for(j=0; j<NHIGGS; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NHIGGS; k++)
           {
           sum+=aux[k]*(B->comp[mf(k,j)]);
           }
        A->comp[mf(i,j)]=sum;
        }
     }
  }


// real part of the trace /NHIGGS
inline double retr_FMatrix(FMatrix const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris;
  double complex tr;

  tr=0.0+0.0*I;
  for(i=0; i<NHIGGS; i++)
     {
     tr+=A->comp[mf(i,i)];
     }
  ris=creal(tr)/(double)NHIGGS;
  return ris;
  }


// imaginary part of the trace /NHIGGS
inline double imtr_FMatrix(FMatrix const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris;
  double complex tr;

  tr=0.0+0.0*I;
  for(i=0; i<NHIGGS; i++)
     {
     tr+=A->comp[mf(i,i)];
     }
  ris=cimag(tr)/(double)NHIGGS;
  return ris;
  }



// l2 norm of the matrix
inline double norm_FMatrix(FMatrix const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris=0.0;

  for(i=0; i<NHIGGS*NHIGGS; i++)
     {
     ris+=cabs(A->comp[i])*cabs(A->comp[i]);
     }

  return sqrt(ris);
  }


#endif
