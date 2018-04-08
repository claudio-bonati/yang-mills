#ifndef U1_C
#define U1_C

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/aligncheck.h"
#include"../include/endianness.h"
#include"../include/macro.h"
#include"../include/random.h"
#include"../include/tens_prod.h"
#include"../include/u1.h"

void init_U1(U1 * restrict A, double complex vec)
  {
  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
    #endif

    #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
    A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
    #else
      #ifdef __INTEL_COMPILER
       __assume_aligned(A, DOUBLE_ALIGN);
      #endif
    #endif
  #endif

  A->comp=vec;
  }


// A=1
void one_U1(U1 * restrict A)
  {
  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
    #endif

    #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
    A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
    #else
      #ifdef __INTEL_COMPILER
       __assume_aligned(A, DOUBLE_ALIGN);
      #endif
    #endif
  #endif

  A->comp=1.0;
  }


// A=0
void zero_U1(U1 * restrict A)
  {
  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
    #endif

    #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
    A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
    #else
      #ifdef __INTEL_COMPILER
       __assume_aligned(A, DOUBLE_ALIGN);
      #endif
    #endif
  #endif

  A->comp=0.0;
  }


// A=B
void equal_U1(U1 * restrict A, U1 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp=B->comp;
  }


// A=B^{dag}
void equal_dag_U1(U1 * restrict A, U1 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp = conj(B->comp);
  }


// A+=B
void plus_equal_U1(U1 * restrict A, U1 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp += B->comp;
  }


// A+=B^{dag}
void plus_equal_dag_U1(U1 * restrict A, U1 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp += conj(B->comp);
  }


// A-=B
void minus_equal_U1(U1 * restrict A, U1 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp -= B->comp;
  }


// A-=(r*B)
void minus_equal_times_real_U1(U1 * restrict A, U1 const * restrict B, double r)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp -= (r*B->comp);
  }


// A-=B^{dag}
void minus_equal_dag_U1(U1 * restrict A, U1 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp -= conj(B->comp);
  }


// A=b*B+c*C
void lin_comb_U1(U1 * restrict A,
                  double b, U1 const * restrict B,
                  double c, U1 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)C, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp = b*B->comp + c*C->comp;
  }


// A=b*B^{dag}+c*C
void lin_comb_dag1_U1(U1 * restrict A,
                       double b, U1 const * restrict B,
                       double c, U1 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)C, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp =  b*conj(B->comp) + c*C->comp;
  }


// A=b*B+c*C^{dag}
void lin_comb_dag2_U1(U1 * restrict A,
                       double b, U1 const * restrict B,
                       double c, U1 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)C, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

       A->comp = b*B->comp + c*conj(C->comp);
  }


// A=b*B^{dag}+c*C^{dag}
void lin_comb_dag12_U1(U1 * restrict A,
                        double b, U1 const * restrict B,
                        double c, U1 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)C, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp = b*conj(B->comp) + c*conj(C->comp);
  }


// A*=r
void times_equal_real_U1(U1 * restrict A, double r)
  {
  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
    #endif

    #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
    A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
    #else
      #ifdef __INTEL_COMPILER
       __assume_aligned(A, DOUBLE_ALIGN);
      #endif
    #endif
  #endif

  A->comp*=r;
  }


// A*=B
void times_equal_U1(U1 * restrict A, U1 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp *= B->comp;
  }


// A*=B^{dag}
void times_equal_dag_U1(U1 * restrict A, U1 const * restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp *= conj(B->comp);
  }


// A=B*C
void times_U1(U1 * restrict A, U1 const * restrict B, U1 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)C, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp= B->comp * C->comp;
  }


// A=B^{dag}*C
void times_dag1_U1(U1 * restrict A, U1 const * restrict B, U1 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)C, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp = conj(B->comp)*C->comp;
  }


// A=B*C^{dag}
void times_dag2_U1(U1 * restrict A, U1 const * restrict B, U1 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)C, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp = B->comp*conj(C->comp);
  }


// A=B^{dag}*C^{dag}
void times_dag12_U1(U1 * restrict A, U1 const * restrict B, U1 const * restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)B, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)C, DOUBLE_ALIGN, __FILE__, __LINE__);
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
  #endif

  A->comp = conj(B->comp)*conj(C->comp);
  }


void rand_matrix_U1(U1 * restrict A)
  {
  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
    #endif

    #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
    A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
    #else
      #ifdef __INTEL_COMPILER
       __assume_aligned(A, DOUBLE_ALIGN);
      #endif
    #endif
  #endif

  double p0, p1, p;

  do
    {
    p0=1.0-2.0*casuale();
    p1=1.0-2.0*casuale();

    p=sqrt(p0*p0+p1*p1);
    }
  while(p<MIN_VALUE);

  p0/=p;
  p1/=p;

  A->comp = p0 + p1*I;
  }


// l2 norm of the matrix
double norm_U1(U1 const * restrict A)
  {
  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
    #endif

    #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
    A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
    #else
      #ifdef __INTEL_COMPILER
       __assume_aligned(A, DOUBLE_ALIGN);
      #endif
    #endif
  #endif

  return sqrt(creal(A->comp)*creal(A->comp)+cimag(A->comp)*cimag(A->comp));
  }


// real part of the trace
double retr_U1(U1 const * restrict A)
  {
  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
    #endif

    #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
    A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
    #else
      #ifdef __INTEL_COMPILER
       __assume_aligned(A, DOUBLE_ALIGN);
      #endif
    #endif
  #endif

  return creal(A->comp);
  }


// imaginary part of the trace
double imtr_U1(U1 const * restrict A)
  {
  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
    #endif

    #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
    A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
    #else
      #ifdef __INTEL_COMPILER
       __assume_aligned(A, DOUBLE_ALIGN);
      #endif
    #endif
  #endif

  return cimag(A->comp);
  }


// unitarize the matrix
void unitarize_U1(U1 * restrict A)
  {
  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
    #endif

    #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
    A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
    #else
      #ifdef __INTEL_COMPILER
       __assume_aligned(A, DOUBLE_ALIGN);
      #endif
    #endif
  #endif

  double p;

  p=norm_U1(A);
  A->comp/=p;
  }

// exponential of the antihermitian part (NO TRACELESS!)
void taexp_U1(U1 * restrict A)
  {
  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)A, DOUBLE_ALIGN, __FILE__, __LINE__);
    #endif

    #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
    A  = __builtin_assume_aligned(A, DOUBLE_ALIGN);
    #else
      #ifdef __INTEL_COMPILER
       __assume_aligned(A, DOUBLE_ALIGN);
      #endif
    #endif
  #endif

  double angle, c, s;

  angle=cimag(A->comp);
  c=cos(angle);
  s=sin(angle);

  A->comp=c+I*s;
  }

// print on screen
void print_on_screen_U1(U1 const * const A)
  {
  printf("%.16lf %.16lf\n", creal(A->comp), cimag(A->comp));
  }


// print on file
void print_on_file_U1(FILE *fp, U1 const * const A)
  {
  int err;
  err=fprintf(fp, "%.16lf %.16lf\n", creal(A->comp), cimag(A->comp));
  if(err<0)
    {
    fprintf(stderr, "Problem in writing on a file a U1 matrix (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// print on binary file without changing endiannes
void print_on_binary_file_noswap_U1(FILE *fp, U1 const * const A)
  {
  size_t err=0;
  double re, im;

  re=creal(A->comp);
  im=cimag(A->comp);

  err=fwrite(&re, sizeof(double), 1, fp);
  if(err!=1)
    {
    fprintf(stderr, "Problem in binary writing on a file a U1 matrix (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  err=fwrite(&im, sizeof(double), 1, fp);
  if(err!=1)
    {
    fprintf(stderr, "Problem in binary writing on a file a U1 matrix (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// print on binary file changing endiannes
void print_on_binary_file_swap_U1(FILE *fp, U1 const * const A)
  {
  double tmp;
  size_t err=0;

  tmp=creal(A->comp);
  SwapBytesDouble(&tmp);
  err+=fwrite(&(tmp), sizeof(double), 1, fp);

  tmp=cimag(A->comp);
  SwapBytesDouble(&tmp);
  err+=fwrite(&(tmp), sizeof(double), 1, fp);

  if(err!=2)
    {
    fprintf(stderr, "Problem in binary writing on a file a U1 matrix (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// print on binary file in big endian format
void print_on_binary_file_bigen_U1(FILE *fp, const U1 * const A)
  {
  if(endian()==0) // little endian machine
    {
    print_on_binary_file_swap_U1(fp, A);
    }
  else
    {
    print_on_binary_file_noswap_U1(fp, A);
    }
  }


// read from file
void read_from_file_U1(FILE *fp, U1 *A)
  {
  double re, im;
  int err;

  err=fscanf(fp, "%lg %lg", &re, &im);
  if(err!=2)
    {
    fprintf(stderr, "Problems reading U1 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  A->comp=re+im*I;
  }


// read from binary file without changing endiannes
void read_from_binary_file_noswap_U1(FILE *fp, U1 *A)
  {
  double re, im;
  size_t err=0;
  err+=fread(&re, sizeof(double), 1, fp);
  err+=fread(&im, sizeof(double), 1, fp);
  if(err!=2)
    {
    fprintf(stderr, "Problems reading U1 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  A->comp=re+I*im;
  }


// read from binary file changing endiannes
void read_from_binary_file_swap_U1(FILE *fp, U1 *A)
  {
  double re, im;
  size_t err=0;

  err+=fread(&re, sizeof(double), 1, fp);
  err+=fread(&im, sizeof(double), 1, fp);

  if(err!=2)
    {
    fprintf(stderr, "Problems reading U1 matrix from file (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  SwapBytesDouble((void *)&re);
  SwapBytesDouble((void *)&im);

  A->comp=re+I*im;
  }


// read from binary file written in big endian
void read_from_binary_file_bigen_U1(FILE *fp, U1 *A)
  {
  if(endian()==0) // little endian machine
    {
    read_from_binary_file_swap_U1(fp, A);
    }
  else
    {
    read_from_binary_file_noswap_U1(fp, A);
    }
  }


void TensProd_init_U1(TensProd * restrict TP, U1 const * restrict A1, U1 const * restrict A2)
  {
  #ifdef MEMALIGN_MODE
    #ifdef DEBUG
      is_aligned((void *)TP, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)A1, DOUBLE_ALIGN, __FILE__, __LINE__);
      is_aligned((void *)A2, DOUBLE_ALIGN, __FILE__, __LINE__);
    #endif

    #if (defined(__GNUC__) && (GCC_VERSION > 40700) ) || defined(__clang__)
    TP = __builtin_assume_aligned(TP, DOUBLE_ALIGN);
    A1 = __builtin_assume_aligned(A1, DOUBLE_ALIGN);
    A2 = __builtin_assume_aligned(A2, DOUBLE_ALIGN);
    #else
      #ifdef __INTEL_COMPILER
       __assume_aligned(TP, DOUBLE_ALIGN);
       __assume_aligned(A1, DOUBLE_ALIGN);
       __assume_aligned(A2, DOUBLE_ALIGN);
      #endif
    #endif
  #endif

  TP->comp[0][0][0][0]=conj(A1->comp)*A2->comp;
  }



#endif
