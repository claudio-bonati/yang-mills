#ifndef SON_H
#define SON_H

#include<complex.h>
#include<math.h>
#include<stdio.h>

#include"flavour_matrix.h"
#include"macro.h"
#include"tens_prod.h"
#include"tens_prod_adj.h"

typedef struct SoN {
   double comp[NCOLOR*NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
} SoN;
//
//  the element [i][j] can be obtained by matrix.comp[m(i,j)] with m(i,j) defined in macro.h
//

typedef struct SoNAdj {
   #if NCOLOR!=1
     double comp[NCOLOR*(NCOLOR-1)/2*NCOLOR*(NCOLOR-1)/2] __attribute__((aligned(DOUBLE_ALIGN)));
   #else // this will never be used, is defined just to avoid warnings
     double comp[1] __attribute__((aligned(DOUBLE_ALIGN)));
   #endif
} SoNAdj;
//
//  the element [i][j] can be obtained by matrix.comp[madj(i,j)] with madj(i,j) defined in macro.h
//

typedef struct SoNVecs {
   double comp[NCOLOR*NHIGGS] __attribute__((aligned(DOUBLE_ALIGN)));
} SoNVecs;
//
// different flavours are associated to different [NCOLOR] blocks of SoNVecs
//

// ***************** for SoN


// A=1
inline void one_SoN(SoN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=0.0;
     }

  for(i=0; i<NCOLOR; i++)
     {
     A->comp[m(i,i)]=1.0;
     }
  }


// A=0
inline void zero_SoN(SoN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=0.0;
     }
  }


// A=B
inline void equal_SoN(SoN * restrict A, SoN const * const restrict B)
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

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=B->comp[i];
     }
  }


// A=B^{dag}
inline void equal_dag_SoN(SoN * restrict A, SoN const * const restrict B)
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

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=B->comp[m(j,i)];
        }
     }
  }


// A+=B
inline void plus_equal_SoN(SoN * restrict A, SoN const * const restrict B)
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

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]+=B->comp[i];
     }
  }


// A+=B^{dag}
inline void plus_equal_dag_SoN(SoN * restrict A, SoN const * const restrict B)
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

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]+=B->comp[m(j,i)];
        }
     }
  }


// A-=B
inline void minus_equal_SoN(SoN * restrict A, SoN const * const restrict B)
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

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]-=B->comp[i];
     }
  }


// A-=(r*B)
inline void minus_equal_times_real_SoN(SoN * restrict A, SoN const * const restrict B, double r)
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

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]-=(r*B->comp[i]);
     }
  }


// A-=B^{dag}
inline void minus_equal_dag_SoN(SoN * restrict A, SoN const * const restrict B)
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

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]-=B->comp[m(j,i)];
        }
     }
  }


// A=b*B+c*C
inline void lin_comb_SoN(SoN * restrict A,
                         double b, SoN const * const restrict B,
                         double c, SoN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=b*(B->comp[i])+c*(C->comp[i]);
     }
  }


// A=b*B^{dag}+c*C
inline void lin_comb_dag1_SoN(SoN * restrict A,
                              double b, SoN const * const restrict B,
                              double c, SoN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=b*B->comp[m(j,i)]+c*(C->comp[m(i,j)]);
        }
     }
  }


// A=b*B+c*C^{dag}
inline void lin_comb_dag2_SoN(SoN * restrict A,
                              double b, SoN const * const restrict B,
                              double c, SoN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=b*(B->comp[m(i,j)])+c*C->comp[m(j,i)];
        }
     }
  }


// A=b*B^{dag}+c*C^{dag}
inline void lin_comb_dag12_SoN(SoN * restrict A,
                               double b, SoN const * const restrict B,
                               double c, SoN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=b*B->comp[m(j,i)]+c*C->comp[m(j,i)];
        }
     }
  }


// A*=r
inline void times_equal_real_SoN(SoN * restrict A, double r)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]*=r;
     }
  }


// A*=r
inline void times_equal_complex_SoN(SoN * restrict A, double complex r)
  {
  (void) A;
  (void) r; // just to avoid warnings

  fprintf(stderr, "The function times_equal_complex_SoN does not exist! (%s, %d)\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);
  }


// A*=B
inline void times_equal_SoN(SoN * restrict A, SoN const * const restrict B)
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
  double aux[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
  double sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux[j]=A->comp[m(i,j)];
        }

     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=aux[k]*(B->comp[m(k,j)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A*=B^{dag}
inline void times_equal_dag_SoN(SoN * restrict A, SoN const * const restrict B)
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
  double aux[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
  double sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux[j]=A->comp[m(i,j)];
        }

     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=aux[k]*B->comp[m(j,k)];
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B*C
inline void times_SoN(SoN * restrict A,
                      SoN const * const restrict B,
                      SoN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[m(i,k)])*(C->comp[m(k,j)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B^{dag}*C
inline void times_dag1_SoN(SoN * restrict A,
                           SoN const * const restrict B,
                           SoN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=B->comp[m(k,i)]*(C->comp[m(k,j)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B*C^{dag}
inline void times_dag2_SoN(SoN * restrict A,
                           SoN const * const restrict B,
                           SoN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[m(i,k)])*(C->comp[m(j,k)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B^{dag}*C^{dag}
inline void times_dag12_SoN(SoN * restrict A,
                            SoN const * const restrict B,
                            SoN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[m(k,i)])*(C->comp[m(j,k)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=lambda*B with lambda diagonal matrix
inline void diag_matrix_times_SoN(SoN * restrict A,
                                  double const lambda[NCOLOR],
                                  SoN const * const restrict B)
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

  int i,j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)] = lambda[i]*(B->comp[m(i,j)]);
        }
     }
  }


// A=lambda*B^{dag} with lambda diagonal matrix
inline void diag_matrix_times_dag_SoN(SoN * restrict A,
                                      double const lambda[NCOLOR],
                                      SoN const * const restrict B)
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

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)] = lambda[i]*B->comp[m(j,i)];
        }
     }
  }


// SO(N) random matrix
// generated a la Cabibbo Marinari with N(N-1)/2 SO(2) random matrices
void rand_matrix_SoN(SoN *A);


// l2 norm of the matrix
inline double norm_SoN(SoN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double aux, ris;

  ris=0.0;
  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     aux=fabs(A->comp[i]);
     ris+=aux*aux;
     }
  return sqrt(ris);
  }


// real part of the trace /N
inline double retr_SoN(SoN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris;
  double tr;

  tr=0.0;
  for(i=0; i<NCOLOR; i++)
     {
     tr+=A->comp[m(i,i)];
     }
  ris=tr/(double)NCOLOR;
  return ris;
  }


// imaginary part of the trace /N
inline double imtr_SoN(SoN const * const restrict A)
  {
  (void) A; // just to avoid warnings

  return 0.0;
  }


// LU decomposition with partial pivoting
void LU_SoN(SoN const * const A, SoN *rif, int *sign);


// determinant
inline double det_SoN(SoN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  #if NCOLOR==2
    double ris=0.0;

    ris+=(A->comp[m(0,0)])*(A->comp[m(1,1)]);
    ris-=(A->comp[m(0,1)])*(A->comp[m(1,0)]);

    return ris;

  #elif NCOLOR==3
    double ris=0.0;

    ris+=(A->comp[m(0,0)])*(A->comp[m(1,1)])*(A->comp[m(2,2)]);
    ris+=(A->comp[m(1,0)])*(A->comp[m(2,1)])*(A->comp[m(0,2)]);
    ris+=(A->comp[m(2,0)])*(A->comp[m(0,1)])*(A->comp[m(1,2)]);
    ris-=(A->comp[m(2,0)])*(A->comp[m(1,1)])*(A->comp[m(0,2)]);
    ris-=(A->comp[m(1,0)])*(A->comp[m(0,1)])*(A->comp[m(2,2)]);
    ris-=(A->comp[m(0,0)])*(A->comp[m(2,1)])*(A->comp[m(1,2)]);

    return ris;
  #else
    int i;
    double ris;
    SoN lu;

    LU_SoN(A, &lu, &i);

    if(i>0)
      {
      ris=1.0;
      }
    else
     {
     ris=-1.0;
     }

    for(i=0; i<NCOLOR; i++)
       {
       ris*=(lu.comp[m(i,i)]);
       }

    return ris;
  #endif
  }


// gives 0 if the matrix is in SO(N) and 1 otherwise
int scheck_SoN(SoN const * const A);


// sunitarize
void unitarize_SoN(SoN *A);


// takes the traceless antihermitian part
inline void ta_SoN(SoN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  SoN aux, aux1;
  double trace;
  int i;

  equal_SoN(&aux, A);
  equal_dag_SoN(&aux1, A);
  minus_equal_SoN(&aux, &aux1);
  times_equal_real_SoN(&aux, 0.5); // now aux=(A-A^{dag})/2

  trace=aux.comp[m(0,0)];
  for(i=1; i<NCOLOR; i++)
     {
     trace+=aux.comp[m(i,i)];
     }
  trace/=(double)NCOLOR;

  for(i=0; i<NCOLOR; i++)
     {
     aux.comp[m(i,i)]-=trace;
     }

  equal_SoN(A, &aux);
  }


// eponential of the traceless antihermitian part
inline void taexp_SoN(SoN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  SoN aux, uno, ris;

  equal_SoN(&aux, A);
  ta_SoN(&aux);

  one_SoN(&uno);

  // now aux is the traceless antihermitian part of the initial matrix
  // and we use
  // exp(x)=1+x(1+x/2(1+x/3*(1+x/4*(1+x/5*....

  equal_SoN(&ris, &aux);
  times_equal_real_SoN(&ris, 0.2);
  plus_equal_SoN(&ris, &uno);

  times_equal_SoN(&ris, &aux);
  times_equal_real_SoN(&ris, 0.25);
  plus_equal_SoN(&ris, &uno);

  times_equal_SoN(&ris, &aux);
  times_equal_real_SoN(&ris, 0.33333333333333333333);
  plus_equal_SoN(&ris, &uno);

  times_equal_SoN(&ris, &aux);
  times_equal_real_SoN(&ris, 0.5);
  plus_equal_SoN(&ris, &uno);

  times_equal_SoN(&ris, &aux);
  plus_equal_SoN(&ris, &uno);

  unitarize_SoN(&ris);
  equal_SoN(A, &ris);
  }



// return 0 if matrix is traceless antihermitian, 1 otherwise
inline int ta_check_SoN(SoN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  double aux;
  int i, j, ris;

  ris=0;

  aux=0.0;
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux+=(A->comp[m(i,j)]+A->comp[m(j,i)]);
        }
     }
  if(fabs(aux)>MIN_VALUE) ris=1;

  aux=0.0;
  for(i=0; i<NCOLOR; i++)
     {
     aux+=A->comp[m(i,i)];
     }
  if(fabs(aux)>MIN_VALUE) ris=1;

  return ris;
  }


// exponential of a TA matrix
inline void exp_of_ta_SoN(SoN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  // we use
  // exp(x)=1+x(1+x/2(1+x/3*(1+x/4*(1+x/5*....

  #ifdef DEBUG
  if(ta_check_SoN(A)!=0)
    {
    fprintf(stderr, "Trying to exp. a non TA matrix! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  SoN aux, uno;

  one_SoN(&uno);
  equal_SoN(&aux, A); // in aux the initial matrix is stored

  times_equal_real_SoN(A, 1.0/5.0);
  plus_equal_SoN(A, &uno);

  times_equal_SoN(A, &aux);
  times_equal_real_SoN(A, 1.0/4.0);
  plus_equal_SoN(A, &uno);

  times_equal_SoN(A, &aux);
  times_equal_real_SoN(A, 1.0/3.0);
  plus_equal_SoN(A, &uno);

  times_equal_SoN(A, &aux);
  times_equal_real_SoN(A, 1.0/2.0);
  plus_equal_SoN(A, &uno);

  times_equal_SoN(A, &aux);
  plus_equal_SoN(A, &uno);

  unitarize_SoN(A);
  }


// print on screen
void print_on_screen_SoN(SoN const * const A);


// print on file
int print_on_file_SoN(FILE *fp, SoN const * const A);


// print on binary file without changing endiannes
int print_on_binary_file_noswap_SoN(FILE *fp, SoN const * const A);


// print on binary file changing endiannes
int print_on_binary_file_swap_SoN(FILE *fp, SoN const * const A);


// print on binary file in bigendian
int print_on_binary_file_bigen_SoN(FILE *fp, SoN const * const A);


// read from file
int read_from_file_SoN(FILE *fp, SoN *A);


// read from binary file without changing endiannes
int read_from_binary_file_noswap_SoN(FILE *fp, SoN *A);


// read from binary file changing endianness
int read_from_binary_file_swap_SoN(FILE *fp, SoN *A);


// read from binary file written in bigendian
int read_from_binary_file_bigen_SoN(FILE *fp, SoN *A);


// initialize tensor product
inline void TensProd_init_SoN(TensProd * restrict TP, SoN const * const restrict A1, SoN const * const restrict A2)
  {
  #ifdef DEBUG
  if(A1==A2)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(A2->comp), DOUBLE_ALIGN);
  __assume_aligned(&(TP->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k, l;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              TP->comp[i][j][k][l]=(A1->comp[m(i,j)])*A2->comp[m(k,l)];
              }
           }
        }
     }
  }


// ***************** for SoNAdj



// convert the fundamental representation matrix B to the adjoint representation matrix A
inline void fund_to_adj_SoN(SoNAdj * restrict A, SoN const * const restrict B)
  {
  (void) A;
  (void) B;

  fprintf(stderr, "The function fund_to_adj_SoN still has to be written (%s, %d)\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);
  }


// initialize tensor product in the adjoint representation
// using two matrices in the fundamental representation
inline void TensProdAdj_init_SoN(TensProdAdj * restrict TP, SoN const * const restrict A1, SoN const * const restrict A2)
  {
  (void) TP;
  (void) A1;
  (void) A2;

  fprintf(stderr, "The function TensProd_adj_init_SoN still has to be written (%s, %d)\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);
  }


// initialize tensor product in the adjoint representation
// using two matrices in the adjoint representation
inline void TensProdAdj_init_SoNAdj(TensProdAdj * restrict TP, SoNAdj const * const restrict A1, SoNAdj const * const restrict A2)
  {
  (void) TP;
  (void) A1;
  (void) A2;

  fprintf(stderr, "The function TensProd_adj_init_SoNAdj still has to be written (%s, %d)\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);
  }


// A=1
inline void one_SoNAdj(SoNAdj * restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*(NCOLOR-1)/2*NCOLOR*(NCOLOR-1)/2; i++)
     {
     A->comp[i]=0.0;
     }

  for(i=0; i<NCOLOR*(NCOLOR-1)/2; i++)
     {
     A->comp[madj(i,i)]=1.0;
     }
  }


// A=0
inline void zero_SoNAdj(SoNAdj * restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*(NCOLOR-1)/2*NCOLOR*(NCOLOR-1)/2; i++)
     {
     A->comp[i]=0.0;
     }
  }


// A+=B
inline void plus_equal_SoNAdj(SoNAdj * restrict A, SoNAdj const * const restrict B)
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
  for(i=0; i<NCOLOR*(NCOLOR-1)/2*NCOLOR*(NCOLOR-1)/2; i++)
     {
     A->comp[i]+=B->comp[i];
     }
  }


// A*=r
inline void times_equal_real_SoNAdj(SoNAdj * restrict A, double r)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*(NCOLOR-1)/2*NCOLOR*(NCOLOR-1)/2; i++)
     {
     A->comp[i]*=r;
     }
  }


// A*=B
inline void times_equal_SoNAdj(SoNAdj * restrict A, SoNAdj const * const restrict B)
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

  #if NCOLOR!=1  // just to avoid warnings
    int i, j, k;
    double aux[NCOLOR*(NCOLOR-1)/2] __attribute__((aligned(DOUBLE_ALIGN)));
    double sum;

    for(i=0; i<NCOLOR*(NCOLOR-1)/2; i++)
       {
       for(j=0; j<NCOLOR*(NCOLOR-1)/2; j++)
          {
          aux[j]=A->comp[madj(i,j)];
          }

       for(j=0; j<NCOLOR*(NCOLOR-1)/2; j++)
          {
          sum=0.0;
          for(k=0; k<NCOLOR*(NCOLOR-1)/2; k++)
             {
             sum+=aux[k]*(B->comp[madj(k,j)]);
             }
          A->comp[madj(i,j)]=sum;
          }
       }
  #else
    (void) A;
    (void) B;
  #endif
  }


// trace in the adjoint rep.
inline double retr_SoNAdj(SoNAdj * restrict A)
  {
  int i;
  double ris=0.0;

  for(i=0; i<NCOLOR*(NCOLOR-1)/2; i++)
     {
     ris+=A->comp[madj(i,i)];
     }

  #if NCOLOR!=1
    ris/=(NCOLOR*(NCOLOR-1)/2);
  #endif

  return ris;
  }


// ***************** for SoNVecs


// A=1
inline void one_SoNVecs(SoNVecs * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NHIGGS; i++)
     {
     A->comp[i]=1.0;
     }
  }


// A=0
inline void zero_SoNVecs(SoNVecs * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NHIGGS; i++)
     {
     A->comp[i]=0.0;
     }
  }


// A=B
inline void equal_SoNVecs(SoNVecs * restrict A, SoNVecs const * const restrict B)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NHIGGS; i++)
     {
     A->comp[i]=B->comp[i];
     }
  }


// A -> A^{\dag}
inline void conjugate_SoNVecs(SoNVecs * restrict A)
  {
  (void) A; // just to avoid warnings
  }


// A-=B
inline void minus_equal_SoNVecs(SoNVecs * restrict A, SoNVecs const * const restrict B)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NHIGGS; i++)
     {
     A->comp[i]-=B->comp[i];
     }
  }


// A+=B
inline void plus_equal_SoNVecs(SoNVecs * restrict A, SoNVecs const * const restrict B)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NHIGGS; i++)
     {
     A->comp[i]+=B->comp[i];
     }
  }


// *= with real number
inline void times_equal_real_SoNVecs(SoNVecs * restrict A, double r)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NHIGGS; i++)
     {
     A->comp[i]*=r;
     }
  }


// *= with real for a single component
inline void times_equal_real_single_SoNVecs(SoNVecs * restrict A, double r, int j)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=NCOLOR*j; i<NCOLOR*(j+1); i++)
     {
     A->comp[i]*=r;
     }
  }


// *= with complex number for a single component
inline void times_equal_complex_single_SoNVecs(SoNVecs * restrict A, double complex r, int j)
  {
  fprintf(stderr, "times_equal_complex_single_SoNVecs can not be used! (%s, %d)\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);

  (void) A;
  (void) r;
  (void) j; // juist to avoid warnings
  }


// norm
inline double norm_SoNVecs(SoNVecs const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris=0.0;

  for(i=0; i<NCOLOR*NHIGGS; i++)
     {
     ris+=(A->comp[i]*A->comp[i]);
     }

  return sqrt(ris);
  }


// normalize
inline void normalize_SoNVecs(SoNVecs * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  double norm=norm_SoNVecs(A);

  times_equal_real_SoNVecs(A, 1.0/norm);
  }


// random vector (normalized)
void rand_vec_SoNVecs(SoNVecs * restrict A);


// real part of the scalar product re(v_1^{\dag}v_2)
inline double re_scal_prod_SoNVecs(SoNVecs const * const restrict v1, SoNVecs const * const restrict v2)
  {
   #ifdef __INTEL_COMPILER
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v2->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris=0.0;

  for(i=0; i<NCOLOR*NHIGGS; i++)
     {
     ris+=(v1->comp[i]) * v2->comp[i];
     }

  return ris;
  }


// real part of the scalar product re(v_1[a]^{\dag}v_2[b]) with a, b flavour indices
inline double re_scal_prod_single_SoNVecs(SoNVecs const * const restrict v1,
                                          SoNVecs const * const restrict v2,
                                          int a,
                                          int b)
  {
   #ifdef __INTEL_COMPILER
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v2->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris=0.0;

  for(i=0; i<NCOLOR; i++)
     {
     ris+=v1->comp[a*NCOLOR+i] * v2->comp[b*NCOLOR+i];
     }

  return ris;
  }


// the i-th component of v2 is multiplied by "matrix"
// v1=matrix*v2
inline void matrix_times_vector_single_SoNVecs(SoNVecs * restrict v1,
                                               SoN const * const restrict matrix,
                                               SoNVecs const * const restrict v2,
                                               int i)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(matrix->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v2->comp), DOUBLE_ALIGN);
  #endif

  int j, k;

  equal_SoNVecs(v1, v2);

  for(j=0; j<NCOLOR; j++)
     {
     v1->comp[NCOLOR*i+j]=0;
     }

  for(j=0; j<NCOLOR; j++)
     {
     for(k=0; k<NCOLOR; k++)
        {
        v1->comp[NCOLOR*i+j] += matrix->comp[m(j,k)] * v2->comp[NCOLOR*i+k];
        }
     }
  }


// all the components of v2 are multiplied by "matrix"
// v1=matrix*v2
inline void matrix_times_vector_all_SoNVecs(SoNVecs * restrict v1,
                                            SoN const * const restrict matrix,
                                            SoNVecs const * const restrict v2)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(matrix->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v2->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        v1->comp[NCOLOR*i+j]=0;
        }

     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           v1->comp[NCOLOR*i+j] += matrix->comp[m(j,k)] * v2->comp[NCOLOR*i+k];
           }
        }
     }
  }


// tensor product of two vectors
// Re(v1^{\dag} * aux * v2) = ReTr(aux * matrix)
inline void vector_tensor_vector_SoNVecs(SoN * restrict matrix,
                                         SoNVecs const * const restrict v1,
                                         SoNVecs const * const restrict v2)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(matrix->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v2->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;

  zero_SoN(matrix);

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           matrix->comp[m(j,k)]+=(v2->comp[NCOLOR*i+j])*(v1->comp[NCOLOR*i+k]);
           }
        }
     }
  }


// initialize the flavour matrix with a vector
// FM[mf(i,j)]=\sum_{on_gauge}conj(v1[i])v1[j] - delta^{ij}/N
// i, j are the flavour indices
inline void init_FMatrix_SoNVecs(FMatrix * restrict fmatrix, SoNVecs const * const restrict v1)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(fmatrix->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;

  zero_FMatrix(fmatrix);

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<NHIGGS; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           fmatrix->comp[mf(i,j)]+=( (v1->comp[NCOLOR*i+k])*(v1->comp[NCOLOR*j+k]) + 0.0*I );
           }
        }
     }

  for(i=0; i<NHIGGS; i++)
     {
     fmatrix->comp[mf(i,i)]-=( 1.0/(double)NHIGGS + 0.0*I);
     }
  }


// return a double coumplex number to check the fate of U(1) flavour symmetry
inline double complex HiggsU1Obs_SoNVecs(SoNVecs const * const restrict v1)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  #endif

  #if NHIGGS >=NCOLOR
    SoN aux;
    int i, j;

    for(i=0; i<NCOLOR; i++)
       {
       for(j=0; j<NCOLOR; j++)
          {
          aux.comp[m(i,j)]=(v1->comp[NCOLOR*i+j]);
          }
       }

    return det_SoN(&aux) + 0.0*I;
  #else
    (void) v1;
    return 0.0 + 0.0*I;
  #endif
  }


// print on file
int print_on_file_SoNVecs(FILE *fp, SoNVecs const * const A);


// print on binary file without changing endiannes
int print_on_binary_file_noswap_SoNVecs(FILE *fp, SoNVecs const * const A);


// print on binary file changing endiannes
int print_on_binary_file_swap_SoNVecs(FILE *fp, SoNVecs const * const A);


// print on binary file in bigendian
int print_on_binary_file_bigen_SoNVecs(FILE *fp, SoNVecs const * const A);


// read from file
int read_from_file_SoNVecs(FILE *fp, SoNVecs *A);


// read from binary file without changing endiannes
int read_from_binary_file_noswap_SoNVecs(FILE *fp, SoNVecs *A);


// read from binary file changing endianness
int read_from_binary_file_swap_SoNVecs(FILE *fp, SoNVecs *A);


// read from binary file written in bigendian
int read_from_binary_file_bigen_SoNVecs(FILE *fp, SoNVecs *A);



#endif // SON_H
